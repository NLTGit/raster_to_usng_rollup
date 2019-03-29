from __future__ import with_statement
import os
import re
import arcpy
import glob
import uuid
import numpy
import string
import random
import hashlib
import sqlite3
import datetime
from arcpy import env
from arcpy.sa import SetNull

from contextlib import closing


# Globals
# TODO Set as tool/gui inputs
ROOTDIR = r'G:\ssdworkspace\raster_to_usng_rollup\data'
WORKSPACE = r'{}\scratch'.format(ROOTDIR)
SQLITEDB = r'{}\db\rollupdata.sqlite'.format(ROOTDIR)
USNGGRID = r'{}\testdata\usng_1k_withpr.shp'.format(ROOTDIR)
#USNGGRID = r'{}\testdata\test.gdb\usng1'.format(ROOTDIR)
ZONEFIELD = "USNG_1KM"
RASTER = r'{}\testdata\Irma_DG_MO_FEMANHRAP_020618_161215.tif'.format(ROOTDIR)
#RASTER = r'{}\testdata\test.gdb\irma1'.format(ROOTDIR)
DISCARDCELLS = None


class DBC:
    """
    The DBC class controls connections and execution to the sqlite database.
    """
    def __init__(self, db):
        """
        Manage the database location and connection

        :param db: Sqlite database location
        :type db: String
        """
        self.db = db
        self.conn = None
        self.create_db()
        self.connect()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """
        Close the connection to the database
        """
        self.conn.commit()
        self.conn.close()

    def create_db(self):
        """
        Create the sqlite db if it doesn't exist
        """
        if not os.path.exists(self.db):
            print("Database does not exist - Creating")
            db_dir = os.path.dirname(self.db)
            if not os.path.exists(db_dir):
                os.makedirs(db_dir)
            arcpy.gp.CreateSQLiteDB(self.db)

    def connect(self):
        """
        Create the database connection
        """
        try:
            self.conn = sqlite3.connect(self.db)
            self.conn.execute("PRAGMA foreign_keys = 1")

        except sqlite3.Error as e:
            print(e)
            self.conn = None

    def execute(self, sql, params=(), returntype=None):
        """
        Execute sql command and give return

        :param sql: A sql statement
        :type sql: String
        :param params: Arguments to pass to the sql string
        :type params: List
        :param returntype: None, 'rowid', 'all', 'one'
        :type returntype: String or None
        :return: The values from the requested query
        :rtype: None or List
        """

        with closing(self.conn.cursor()) as cursor:
            try:
                cursor.execute(sql, params)
                if not returntype:
                    return
                elif returntype == 'rowid':
                    return [cursor.lastrowid]
                elif returntype == 'all':
                    return cursor.fetchall()
                elif returntype == 'one':
                    return cursor.fetchone()

            except sqlite3.Error as e:
                print(e)
                return

    def executemany(self, sql, params=()):
        """
        Execute sql command and give return

        :param sql: A sql statement
        :type sql: String
        :param params: Arguments to pass to the sql string
        :type params: List
        """
        with closing(self.conn.cursor()) as cursor:
            try:
                cursor.executemany(sql, params)

            except sqlite3.Error as e:
                print(e)
                return


class RasterProcessor:
    """
    Management of the zonal stats database
    """
    def __init__(self, db, poly, zonefield, raster):
        self.db = db
        self.poly = poly
        self.polyname = os.path.basename(poly)
        self.polyhash = self.generate_md5(poly)
        self.zonefield = zonefield
        self.raster = raster
        self.rastername = os.path.basename(raster)
        self.rasterhash = self.generate_md5(raster)
        self.null_raster = None
        self.zonal_stats_data = []
        self.timestamp = self.timestamp_from_rastername()

    def process_raster(self):
        if self.precheck_db():
            print("Warning: Identical data found in database")
            print("Discontinuing processing")
            return

        self.get_zonal_stats_np_array(self.raster)
        self.update_db()

    def precheck_db(self):
        """
        Creates the schema and determines if analysis settings have been used

        :return: Whether the analysis settings have been used before
        :rtype: Bool
        """
        self.create_sqlite_schema()
        return self.analysis_has_been_run()

    def update_db(self):
        rasterid = self.load_raster_table()
        self.load_poly_table(rasterid)
        self.load_zonal_table(rasterid)
        self.load_footprint_table(rasterid)

    def create_sqlite_schema(self):
        """
        Create the schema for the db
        """

        with DBC(self.db) as dbc:
            dbc.execute("""
                CREATE TABLE IF NOT EXISTS raster (
                  id INTEGER PRIMARY KEY,
                  name TEXT NOT NULL,
                  hash TEXT NOT NULL,
                  referencedate TEXT);""")

            dbc.execute("""CREATE TABLE IF NOT EXISTS poly (
                  raster_id INTERGER NOT NULL,
                  name TEXT NOT NULL,
                  hash TEXT NOT NULL,
                  zone TEXT,
                  CONSTRAINT fk_raster
                    FOREIGN KEY (raster_id)
                    REFERENCES raster(id)
                    ON DELETE CASCADE);""")

            dbc.execute("""CREATE TABLE IF NOT EXISTS raster_footprint (
                  raster_id INTEGER NOT NULL,
                  the_geom TEXT,
                  CONSTRAINT fk_raster
                    FOREIGN KEY (raster_id)
                    REFERENCES raster(id)
                    ON DELETE CASCADE);""")

            dbc.execute("""CREATE TABLE IF NOT EXISTS zonal_stats (
                  raster_id INTEGER NOT NULL,
                  feature_id TEXT,
                  count INTEGER,
                  min REAL,
                  max REAL,
                  mean REAL,
                  median REAL,
                  std REAL,
                  CONSTRAINT fk_raster
                    FOREIGN KEY (raster_id)
                    REFERENCES raster(id)
                    ON DELETE CASCADE);""")

    def load_poly_table(self, rasterid):
        """
        Populate the polygon table used to process each raster

        :param rasterid: Row id corresponding to the raster table
        :type rasterid: Int
        """
        with DBC(self.db) as dbc:
            insert_sql = """
                INSERT into poly(raster_id, name, hash, zone)
                VALUES (?,?,?,?)
            """
            dbc.execute(insert_sql, (rasterid, self.polyname, self.polyhash, self.zonefield))

    def load_zonal_table(self, rasterid):
        """
        Populate the zonal stats table from the object data

        :param rasterid: Row id corresponding to the raster table
        :type rasterid: Int
        """

        with DBC(self.db) as dbc:
            insert_sql = """
                INSERT into zonal_stats(raster_id, feature_id, count, min, max, mean, median, std)
                VALUES (?,?,?,?,?,?,?,?)
            """
            process_list = [[rasterid] + list(row) for row in self.zonal_stats_data]
            dbc.executemany(insert_sql, (process_list))

    def load_raster_table(self):
        """
        Load the raster name and hash into the db, give the sqlite rowid

        :return: rasterid
        :rtype: Integer
        """
        with DBC(self.db) as dbc:
            insert_sql = """
                INSERT INTO raster(name, hash, referencedate) 
                VALUES (?,?,?)
                """
            rowid = dbc.execute(insert_sql, (self.rastername, self.rasterhash, self.timestamp), 'rowid')
            select_sql = """
                SELECT id FROM raster WHERE rowid = ?
                """
            rasterid = dbc.execute(select_sql, (rowid), 'one')

            return rasterid[0]

    def load_footprint_table(self, rasterid):
        """
        Create and load raster WKT extents into the db by unique id

        :param rasterid: Unique raster id
        :type rasterid: Integer
        """

        # get the raster extent and spatial reference
        arcpy_extent = str(arcpy.Raster(self.raster).extent)
        arcpy_in_sref = arcpy.Describe(self.raster).spatialReference
        XMin, YMin, XMax, YMax = arcpy_extent.split()[:4]

        # Use the extent to create a polygon, project to wgs84
        point_array = arcpy.Array([
                             arcpy.Point(XMin, YMin),
                             arcpy.Point(XMin, YMax),
                             arcpy.Point(XMax, YMax),
                             arcpy.Point(XMax, YMin),
                             arcpy.Point(XMin, YMin)
                             ])
        arcpy_out_sref = arcpy.SpatialReference('Geographic Coordinate Systems/World/WGS 1984')
        polygon = arcpy.Polygon(point_array, arcpy_in_sref).projectAs(arcpy_out_sref)

        # update database
        with DBC(self.db) as dbc:
            insert_sql = """
                INSERT INTO raster_footprint(raster_id, the_geom) 
                VALUES (?,?)
                """
            dbc.execute(insert_sql, (rasterid, polygon.WKT))

    def truncate_db(self):
        with DBC(self.db) as dbc:
            delete_sql = """
                DELETE FROM raster;
            """
            dbc.execute(delete_sql)

    def clear_db(self):
        with DBC(self.db) as dbc:
            drop_sql = """DROP TABLE raster;"""
            dbc.execute(drop_sql)

            drop_sql = """DROP TABLE poly;"""
            dbc.execute(drop_sql)

            drop_sql = """DROP TABLE raster_footprint;"""
            dbc.execute(drop_sql)

            drop_sql = """DROP TABLE zonal_stats;"""
            dbc.execute(drop_sql)

    def analysis_has_been_run(self):
        """
        Check that the analysis hasn't been run before with the same config

        :return: Whether the raster and polygon are in the db
        :rtype: Bool
        """
        with DBC(self.db) as dbc:
            # check whether the hash values and column zone key exist in the DB
            select_sql = """
                SELECT 1
                FROM raster r
                LEFT JOIN poly p 
                ON r.id = p.raster_id
                WHERE r.hash = ? AND
                      p.hash = ? AND
                      p.zone = ?;
            """
            results = dbc.execute(select_sql,
                                 (self.rasterhash,
                                  self.polyhash,
                                  self.zonefield),
                                  'one')
            if results:
                return True
            else:
                return False

    def generate_md5(self, file_name):
        """
        Generate the md5sum of an input file without reading it all to memory

        :param file_name: Input file to be hashed, must be readable
        :type file_name: String
        :return: Hash value for file at location
        :rtype: String
        """
        hash_md5 = hashlib.md5()

        # if in a gdb, handle the cases
        if ".gdb" in file_name:
            desc=arcpy.Describe(file_name)

            # raster gdb, hash the file size and band standard dev
            if desc.dataType == "RasterDataset":
                rast=arcpy.Raster(file_name)
                stdev = arcpy.GetRasterProperties_management(rast, "STD")
                hash_md5.update(str(rast.uncompressedSize)+str(stdev))

            # shp gdb, hash the extent and feature count
            if desc.dataType == "FeatureClass":
                fcount = arcpy.GetCount_management(file_name).getOutput(0)
                desc = arcpy.Describe(file_name)
                hash_md5.update(str(fcount)+str(desc.extent))

            else:
                print("ERROR: Unrecognized gdb input type {}, not feature or raster".format(desc.dataType))
                print("File: {}".format(file_name))

        # if the file is a shapefile on the fs, hash the shp and dbf, combine
        elif file_name[-4:] == '.shp' and os.path.exists(file_name[:-4] + ".dbf"):
            # hash the shp
            hash_shp = hashlib.md5()
            with open(file_name, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_shp.update(chunk)

            # hash the dbf
            hash_dbf = hashlib.md5()
            dbf_name = file_name[:-4] + ".dbf"
            with open(dbf_name, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_dbf.update(chunk)

            # hash both files together
            hash_md5.update(hash_shp.hexdigest() + hash_shp.hexdigest())

        # else if the file is on the fs, hash it
        else:
            with open(file_name, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)

        return hash_md5.hexdigest()

    def process_null_raster(self):
        """
        Remove cells from a raster by expression. 'Value < 0' for example

        :param input_raster_path: Path of the raster to be processed
        :type input_raster_path: String
        :param discard_expression: Expression of values to be removed
        :type discard_expression: String
        :return: Path to location of corrected raster
        :rtype: String
        """
        # Discard all raster cells where the value is <= 0
        raster_desc = arcpy.Describe(self.raster)
        discard_condition = "VALUE < 0"
        null_raster = SetNull(self.raster, self.raster, discard_condition)

        # Save the output
        raster_path, raster_file = os.path.split(RASTER)
        raster_name, raster_ext = os.path.splitext(raster_file)
        if raster_ext:
            null_raster_name = "{}_null{}".format(raster_name, raster_ext)
        else:
            null_raster_name = "{}_null".format(raster_name)

        null_raster_path = os.path.join(raster_path, null_raster_name)
        null_raster.save(null_raster_path)
        self.null_raster = null_raster_path

    def generate_random_file(self, ext):
        """

        create a 10 digit random filename for a given ext with a
        non numeric first

        :param ext: Extension
        :type ext: String
        :return: A random filename
        :rtype: String
        """
        # if begins with . discard
        if ext[0] == ".":
            ext = ext[1:]

        random_letter = random.choice(string.ascii_lowercase)
        random_uuid = str(uuid.uuid4().hex)[:9]
        return "{}{}.{}".format(random_letter, random_uuid, ext)

    def get_zonal_stats_np_array(self, raster):
        """
        Generate the zonal stats of an input raster file

        :param raster: The location of a raster on disk
        :type raster: String
        """
        arcpyraster = arcpy.Raster(raster)
        raster_desc = arcpy.Describe(arcpyraster)

        # necessary raster metadata
        oid_fieldname = arcpy.Describe(self.poly).OIDFieldName
        lowerLeft = arcpy.Point(arcpyraster.extent.XMin, raster_desc.extent.YMin)
        cellSize = arcpyraster.meanCellWidth
        raster_nodataval = raster_desc.noDataValue

        # environment settings (some of these may not apply, or cellsize overrides(
        arcpy.env.outputCoordinateSystem = raster_desc.spatialReference
        arcpy.env.snapRaster = arcpyraster
        arcpy.env.extent = arcpyraster.extent

        # Process: Polygon to Raster
        temp_raster = self.generate_random_file("tif")
        arcpy.PolygonToRaster_conversion(self.poly,
                                         oid_fieldname,
                                         temp_raster,
                                         "CELL_CENTER",
                                         "NONE",
                                         cellSize)

        # Clip the raster a second time to verify that the extents match
        # Some cases exist where the cellsize gives an extra data column
        # This is a workaround
        temp_raster2 = self.generate_random_file("tif")
        extent_text = " ".join(str(arcpyraster.extent).split()[:4])
        arcpy.Clip_management(temp_raster,
                              extent_text,
                              temp_raster2,
                              arcpyraster,
                              maintain_clipping_extent="NO_MAINTAIN_EXTENT")

        #reset env
        arcpy.env.outputCoordinateSystem = None
        arcpy.env.snapRaster = None
        arcpy.env.extent = None

        # load the rasters via numpy
        np_poly_raster = arcpy.RasterToNumPyArray(temp_raster2,lowerLeft, nodata_to_value=-1)
        np_raster = arcpy.RasterToNumPyArray(arcpyraster,lowerLeft, nodata_to_value=raster_nodataval)

        # cleanup temp filesystem rasters after conversion
        self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster)+"*",glob_opt=True)
        self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster2)+"*",glob_opt=True)

        if np_raster.shape != np_poly_raster.shape:
            print("ERROR: Rasters for stats calc have differing shape")

        # get unique raster values, remove nodata value
        poly_raster_vals = numpy.unique(np_poly_raster).tolist()
        if -1 in poly_raster_vals: poly_raster_vals.remove(-1)

        # join back the zone name string to the OID, keep as lookup
        oid_zone_lookup = {}
        with arcpy.da.SearchCursor(self.poly, [oid_fieldname, self.zonefield]) as cursor:
            for row in cursor:
                # if OID occurrs in raster from polygon, add to lookup table
                if row[0] in poly_raster_vals:
                    oid_zone_lookup[row[0]] = row[1]


        # create statistics on the input raster by unique values
        for oid in oid_zone_lookup.keys():
            # filter by unique raster value (e.x. 2)
            statsdata_wnodata = np_raster[(np_poly_raster == oid)]
            # filter again and remove nodata values from the raster
            statsdata = statsdata_wnodata[statsdata_wnodata != raster_nodataval]
            if statsdata.size != 0:
                self.zonal_stats_data.append([oid_zone_lookup[oid],
                                              statsdata.size,
                                              numpy.min(statsdata).item(),
                                              numpy.max(statsdata).item(),
                                              numpy.mean(statsdata).item(),
                                              numpy.median(statsdata).item(),
                                              numpy.std(statsdata).item()])

    def arcpy_rm_file(self, input_file, glob_opt=False):
        """
        Removes a file or glob pattern using the arcgis engine.

        This prevents conflicts where arc has a file lock while
        the process cleans up.

        :param input_file: File path or pattern
        :type input_file: String
        :param glob_opt: Flag to notify function of pattern
        :type glob_opt: Bool
        :return: None
        :rtype: None
        """
        try:
            if glob_opt:
                for match in glob.glob(input_file):
                    arcpy.Delete_management(match)
            else:
                arcpy.Delete_management(input_file)
        except OSError:
            pass

    def timestamp_from_rastername(self):
        """
        Generate a timestamp in expected format based on input filename

        :return: Timestamp
        :rtype: String
        """

        pattern_2digityear = "[0-9]{6}_[0-9]{6}" # 020618_161215
        pattern_4digityear = "[0-9]{8}_[0-9]{6}" # 02062018_161215

        result = re.search(pattern_2digityear, self.rastername)
        if not result:
            result = re.search(pattern_4digityear, self.rastername)
        if not result:
            return None

        timestamp_text = result.group()
        ymd, hms = timestamp_text.split("_")
        if len(ymd) == 6:
            dd = datetime.datetime.strptime(timestamp_text, '%m%d%y_%H%M%S')
            return dd.strftime('%Y-%m-%d %H:%M:%S')
        elif len(ymd) == 8:
            dd = datetime.datetime.strptime(timestamp_text, '%m%d%Y_%H%M%S')
            return dd.strftime('%Y-%m-%d %H:%M:%S')

if __name__ == "__main__":

    # Config
    env.overwriteOutput = 1
    env.workspace = WORKSPACE

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    proc = RasterProcessor(SQLITEDB, USNGGRID, ZONEFIELD, RASTER)
    proc.clear_db()
    #proc.truncate_db()
    proc.process_raster()
