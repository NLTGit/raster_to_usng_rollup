from __future__ import with_statement
import os
import re
import sys
import arcpy
import glob
import uuid
import numpy
import string
import random
import hashlib
import sqlite3
import datetime
from contextlib import closing


class DBC:
    """
    The DBC class controls connections and execution to the sqlite database.
    """
    def __init__(self, db, messages):
        """
        Manage the database location and connection

        :param db: Sqlite database location
        :type db: String
        """
        self.db = db
        self.messages = messages
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
            self.messages.addMessage("Database does not exist - Creating")
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
            self.messages.addError(e)
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
                self.messages.addError(e)
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
                self.messages.addError(e)
                return


class RasterProcessor:
    """
    Management of the zonal stats database
    """
    def __init__(self, db, inpoly, outpoly, outfootprints, zonefield, raster, start_timestamp=None, end_timestamp=None, messages=None):
        self.db = db
        self.poly = inpoly
        self.messages = messages
        self.polyname = os.path.basename(inpoly)
        self.polyhash = self.generate_md5(inpoly)
        self.outpoly = outpoly
        self.zonefield = zonefield
        self.raster = raster
        self.rastername = os.path.basename(raster)
        self.rasterhash = self.generate_md5(raster)
        self.null_raster = None
        self.zonal_stats_data = []
        self.timestamp = self.timestamp_from_text(self.rastername)
        self.start_timestamp = self.validate_timestamp(start_timestamp)
        self.end_timestamp = self.validate_timestamp(end_timestamp)
        self.outfootprints = outfootprints


    def process_raster(self):
        self.messages.addWarningMessage("processing raster: {}".format(self.rastername))
        if self.precheck_db():
            self.messages.addWarningMessage("Warning: Identical data found in database")
            self.messages.addWarningMessage("Shortcutting processing for this raster")
            return

        self.get_zonal_stats_np_array(self.raster)
        self.update_db()

    def create_output(self):
        """
        Create the output file and append all relevant data
        """
        self.messages.addMessage("creating output files")
        self.create_output_file()
        self.create_output_footprints()
        for rasterid in self.get_current_poly_rasterids():
            self.join_zonalstats_by_rasterid(rasterid)
            self.join_footprints_by_rasterid(rasterid)

    def precheck_db(self):
        """
        Creates the schema and determines if analysis settings have been used

        :return: Whether the analysis settings have been used before
        :rtype: Bool
        """
        self.messages.addMessage("checking if previous analysis settings have been used before")
        self.create_sqlite_schema()
        return self.analysis_has_been_run()

    def update_db(self):
        self.messages.addMessage("updating local sqlite database with poly-raster data")
        rasterid = self.load_raster_table()
        self.load_poly_table(rasterid)
        self.load_zonal_table(rasterid)
        self.load_footprint_table(rasterid)

    def create_sqlite_schema(self):
        """
        Create the schema for the db
        """
        self.messages.addMessage("creating sqlite table schema")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("populating the polygon table data")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("populating the zonalstatistics table")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("populating the raster table")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("populating the raster footprints table")
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
        with DBC(self.db, self.messages) as dbc:
            insert_sql = """
                INSERT INTO raster_footprint(raster_id, the_geom) 
                VALUES (?,?)
                """
            dbc.execute(insert_sql, (rasterid, polygon.WKT))

    def get_current_poly_rasterids(self):
        """
        Get the raster ids that are linked with this specific polygon. Order results by date timestamp, earliest first. Filter the results by the start and end dates if exists.

        :return: Raster ids in ascending order by timestamp
        :rtype: List
        """
        self.messages.addMessage("getting all relevant rasters in the db analyzed by this polygon")
        with DBC(self.db, self.messages) as dbc:
            # if the timestamp filters are active, filter the data appearing in the output file
            if self.start_timestamp and self.end_timestamp:
                select_sql = """
                  SELECT id
                  FROM raster r
                  WHERE r.id IN (
                    SELECT raster_id
                    FROM poly p
                    WHERE p.name = ? AND
                          p.hash = ? )
                    AND r.referencedate BETWEEN ? AND ?
                  ORDER BY r.referencedate ASC;
                """
                results = dbc.execute(select_sql, (self.polyname, self.polyhash, self.start_timestamp, self.end_timestamp), 'all')

            # otherwise use all data paired with the polygon, ignoring time
            else:
                select_sql = """
                  SELECT id
                  FROM raster r
                  WHERE r.id IN (
                    SELECT raster_id
                    FROM poly p
                    WHERE p.name = ? AND
                          p.hash = ? )
                  ORDER BY r.referencedate ASC;
                """
                results = dbc.execute(select_sql, (self.polyname, self.polyhash), 'all')

            return [r[0] for r in results]

    def truncate_db(self):
        self.messages.addWarningMessage("truncating database")
        with DBC(self.db, self.messages) as dbc:
            delete_sql = """
                DELETE FROM raster;
            """
            dbc.execute(delete_sql)

    def clear_db(self):
        self.messages.addWarningMessage("clearing database tables")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("determine if identical analysis has already been run")
        with DBC(self.db, self.messages) as dbc:
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
        self.messages.addMessage("generating md5sum for {}".format(file_name))
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
                self.messages.addError("Unrecognized gdb input type {}, not feature or raster".format(desc.dataType))
                self.messages.addError("File: {}".format(file_name))

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

    def generate_random_file(self, ext):
        """

        create a 10 digit random filename for a given ext with a
        non numeric first

        :param ext: Extension
        :type ext: String
        :return: A random filename
        :rtype: String
        """
        self.messages.addMessage("generating random file with extension: {}".format(ext))

        random_letter = random.choice(string.ascii_lowercase)
        random_uuid = str(uuid.uuid4().hex)[:9]
        if ext == "":
            return "{}{}".format(random_letter, random_uuid)
        else:
            # if begins with . discard
            if ext[0] == ".":
                ext = ext[1:]
            return "{}{}.{}".format(random_letter, random_uuid, ext)

    def get_zonal_stats_np_array(self, raster):
        """
        Generate the zonal stats of an input raster file

        :param raster: The location of a raster on disk
        :type raster: String
        """
        self.messages.addMessage("generating zonal statistics")
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

        # If workspace is a folder, make a tif, else put in gdb
        if arcpy.env.workspace[-4:].lower() == '.gdb':
            temp_raster = self.generate_random_file("")
            temp_raster2 = self.generate_random_file("")
        else:
            temp_raster = self.generate_random_file("tif")
            temp_raster2 = self.generate_random_file("tif")

        # Process polygon to raster
        arcpy.PolygonToRaster_conversion(self.poly,
                                         oid_fieldname,
                                         temp_raster,
                                         "CELL_CENTER",
                                         "NONE",
                                         cellSize)

        # Clip the raster a second time to verify that the extents match
        # Some cases exist where the cellsize gives an extra data column
        # This is a workaround
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
        np_poly_raster = arcpy.RasterToNumPyArray(temp_raster2, lowerLeft, nodata_to_value=-1)
        np_raster = arcpy.RasterToNumPyArray(arcpyraster, lowerLeft, nodata_to_value=raster_nodataval)

        # there's a bug in arcgis 10 related to background geoprocessing where snap raster adds 1 extra row
        # if there's an extra row in the numpy array, remove it.
        # https://gis.stackexchange.com/questions/34085/aligning-two-non-coincident-equi-resolution-raster-grids-in-arcgis-desktop/34172#34172
        if np_poly_raster.shape[0]-1 == np_raster.shape[0]:
            np_poly_raster = np_poly_raster[:-1, :]

        if np_raster.shape != np_poly_raster.shape:
            self.messages.addWarningMessage("Rasters for stats calc have differing shape")
            self.messages.addMessage("np_raster.shape: {}".format(np_raster.shape))
            self.messages.addMessage("np_poly_raster.shape: {}".format(np_poly_raster.shape))
            arcpy.AddError("Mismatched Raster Dims")
            sys.exit(1)

        # cleanup temp filesystem rasters after conversion
        if arcpy.env.workspace[-4:].lower() == '.gdb':
            self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster))
            self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster2))
        else:
            self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster)+"*",glob_opt=True)
            self.arcpy_rm_file(os.path.join(arcpy.env.workspace,temp_raster2)+"*",glob_opt=True)

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
        self.messages.addMessage("removing input file {}".format(input_file))
        try:
            if glob_opt:
                for match in glob.glob(input_file):
                    arcpy.Delete_management(match)
            else:
                arcpy.Delete_management(input_file)
        except OSError:
            pass

    def validate_timestamp(self, timestamp):
        """
        Validate that a timestamp is y-m-d h:m:s or that it follows the expected
        input datestamp similar to what's seen on the raster file names

        Input can be: y-m-d h:m:s, mmddyy_hhmmss, or mmddyyyy_hhmmss

        :param timestamp: An input timestamp
        :type timestamp: String
        :return: Formatted timestamp or None
        :rtype: String
        """
        # shortcut no input on the timestamp (default handling for class)
        if not timestamp:
            return None
        self.messages.addMessage("validating timestamp {}".format(timestamp))
        # check if input timestamp follows raster name format
        pattern_2digityear = "[0-9]{6}_[0-9]{6}" # 020618_161215
        pattern_4digityear = "[0-9]{8}_[0-9]{6}" # 02062018_161215
        if not re.search(pattern_2digityear, timestamp) \
                and not re.search(pattern_4digityear, timestamp):

            # check if input timestamp is y-m-d h:m:s
            try:
                datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
                return timestamp
            except ValueError:
                self.messages.addError("Incorrect date format should be y-m-d h:m:s, mmddyy_hhmmss, or mmddyyyy_hhmmss")
                return None
        else:
            # convert to y-m-d h:m:s
            return self.timestamp_from_text(timestamp)

    def timestamp_from_text(self, intext=None):
        """
        Generate a timestamp in expected format based on input filename

        :return: Timestamp
        :rtype: String
        """

        self.messages.addMessage("generating timestamp from text: {}".format(intext))

        pattern_2digityear = "[0-9]{6}_[0-9]{6}" # 020618_161215
        pattern_4digityear = "[0-9]{8}_[0-9]{6}" # 02062018_161215

        result = re.search(pattern_2digityear, intext)
        if not result:
            result = re.search(pattern_4digityear, intext)
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

    def create_output_file(self):
        """
        Crete the output shapefile to be joined with statistics
        """
        self.messages.addMessage("creating output polygon file")
        arcpy.CopyFeatures_management(self.poly, self.outpoly)
        try:
            arcpy.RemoveIndex_management(self.poly, ["id_idx"])
        except:
            pass
        arcpy.AddIndex_management(self.outpoly, self.zonefield, "id_idx")

    def create_output_footprints(self):
        """
        Create the footprints output table showing the raster coverage and collection dates
        """
        self.messages.addMessage("creating output footprint file")
        out_path = os.path.dirname(self.outfootprints)
        out_name = os.path.basename(self.outfootprints)
        out_sref = arcpy.SpatialReference('Geographic Coordinate Systems/World/WGS 1984')
        arcpy.CreateFeatureclass_management(out_path, out_name, 'POLYGON', spatial_reference=out_sref)
        arcpy.AddField_management(self.outfootprints, "name", "TEXT", 255)
        arcpy.AddField_management(self.outfootprints, "date", "DATE")
        arcpy.DeleteField_management (self.outfootprints, "id")

    def join_zonalstats_by_rasterid(self, rasterid):
        """
        Joins the zonal stats data to a copy of the input polygon.
        Attribute table is wide so columns are added as needed

        :param rasterid: The rasterid of a dataset to be added
        :type rasterid: Integer
        """
        self.messages.addMessage("joining the zonal statistics to the output polygon")
        # keep running list of fields to be updated
        update_fields = [self.zonefield]

        # create necessary fields
        arcpy.AddField_management(self.outpoly, "count_{}".format(rasterid), "LONG")
        update_fields.append("count_{}".format(rasterid))
        for fieldname in ["min","max","mean","median","std"]:
            fieldname_wid = "{}_{}".format(fieldname, rasterid)
            update_fields.append(fieldname_wid)
            arcpy.AddField_management(self.outpoly, fieldname_wid, "FLOAT")

        with DBC(self.db, self.messages) as dbc:
            select_sql = """
                SELECT feature_id, count, min, max, mean, median, std
                FROM zonal_stats z
                WHERE z.raster_id = ?
            """
            results = dbc.execute(select_sql, (rasterid,), 'all')
        results_dict = {r[0]: r[1:] for r in results}

        def quote_str(instr):
            return "'"+instr+"'"
        arcpy_sql = '"{0}" IN ({1})'.format(self.zonefield,
                                            ', '.join(map(quote_str, results_dict.keys())) or 'NULL')

        # select all features with rasterstat data and update using dictionary
        # FIELDS: id_field + new count, min, max, mean, median, std
        with arcpy.da.UpdateCursor(self.outpoly, update_fields, arcpy_sql) as cursor:
            for row in cursor:
                feature_id = row[0] # id field
                row[1] = results_dict[feature_id][0] # count
                row[2] = results_dict[feature_id][1] # min
                row[3] = results_dict[feature_id][2] # max
                row[4] = results_dict[feature_id][3] # mean
                row[5] = results_dict[feature_id][4] # median
                row[6] = results_dict[feature_id][5] # std
                cursor.updateRow(row)

    def join_footprints_by_rasterid(self, rasterid):
        """
        Join the polygon, name, and date fields to the output raster polygon shapefile by rasterid

        :param rasterid: the rasterid of a dataset to be appended
        :type rasterid: Integer
        """
        self.messages.addMessage("join footprints to the output polygon")
        with DBC(self.db, self.messages) as dbc:
            select_sql = """
                SELECT r.the_geom, l.name, l.referencedate
                FROM raster l
                LEFT JOIN raster_footprint r on l.id = r.raster_id
                WHERE l.id = ?
            """
            results = dbc.execute(select_sql, (rasterid,), 'one')
            with arcpy.da.InsertCursor(self.outfootprints, ['SHAPE@WKT', 'name', 'date']) as cursor:
                cursor.insertRow(results)


class Toolbox(object):
    def __init__(self):
        """Define the toolbox"""
        self.label = "usng_rollup"
        self.alias = "USNG Rollup"

        # List of tool classes associated with this toolbox
        self.tools = [USNGRollup]


class USNGRollup(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Process Data"
        self.description = "This tool generates raster statistics for a set of input rasters by an input polygon"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # work directory
        param0 = arcpy.Parameter(
            displayName="Workspace",
            name="WORKSPACE",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        # In the tool's dialog box, the first parameter will show
        #  the workspace environment's value (if set)
        param0.defaultEnvironmentName = "workspace"

        param1 = arcpy.Parameter(
            displayName="SQLite database",
            name="SQLITEDB",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param1.value = "db.sqlite"

        # Polygon feature to perform rollup
        param2 = arcpy.Parameter(
            displayName="Input polygon (e.x. USNG shp)",
            name="USNGGRID",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        param2.filter.list = ["Polygon"]

        # Name field to be stored in the DB
        param3 = arcpy.Parameter(
            displayName="Unique field name (e.x. USNG name column)",
            name="ZONEFIELD",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param3.parameterDependencies = [param2.name]

        # Raster lists
        param4 = arcpy.Parameter(
            displayName="Input rasters",
            name="RASTER_LIST",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input",
            multiValue=True)

        param5 = arcpy.Parameter(
            displayName="Start timestamp filter (optional) - mmddyy_hhmmss, mmddyyyy_hhmmss, or yyyy-mm-dd hh:mm:ss",
            name="STARTTIMESTAMP",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        param6 = arcpy.Parameter(
            displayName="End timestamp filter (optional) - mmddyy_hhmmss, mmddyyyy_hhmmss, or yyyy-mm-dd hh:mm:ss",
            name="ENDTIMESTAMP",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        # Output statistics polygon - a copy of the USNG file with attached statistics
        param7 = arcpy.Parameter(
            displayName="Output statistics shp or feature class",
            name="OUTPUTGRID",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")
        param7.parameterDependencies = [param2.name]
        param7.value = "statistics.shp"

        # Output footprints polygon
        param8 = arcpy.Parameter(
            displayName="Output raster footprints shp or feature class",
            name="OUTPUTFOOTPRINTS",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")
        param8.parameterDependencies = [param2.name]
        param8.value = "footprints.shp"

        params = [param0, param1, param2, param3, param4, param5, param6, param7, param8]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        WORKSPACE = parameters[0].valueAsText
        RASTER_LIST = parameters[4].valueAsText.split(";")
        SQLITEDB = parameters[1].valueAsText
        USNGGRID = parameters[2].valueAsText
        OUTPUTGRID = parameters[7].valueAsText
        OUTPUTFOOTPRINTS = parameters[8].valueAsText
        ZONEFIELD = parameters[3].valueAsText
        STARTTIMESTAMP = parameters[5].valueAsText
        ENDTIMESTAMP = parameters[6].valueAsText

        # == SQLITE DB== #
        # if the extension isn't .sqlite, add it
        if os.path.splitext(SQLITEDB)[1] != '.sqlite':
            _, file_extension = os.path.splitext(SQLITEDB)
            SQLITEDB = SQLITEDB[:len(SQLITEDB)-len(file_extension)]+".sqlite"
            messages.addWarningMessage("replaced database extension in filename")
            messages.addMessage(SQLITEDB)

        # if the full filepath of the database hasn't been specified
        if not os.path.dirname(SQLITEDB):

            # put it in the workspace parent dir if geodatabase
            if os.path.splitext(WORKSPACE)[1].lower() == ".gdb":
                messages.addWarningMessage("Workspace is not a folder - placing sqlite db in workspace parent directory")
                parent_dir = os.path.dirname(WORKSPACE)
                SQLITEDB = os.path.join(parent_dir, SQLITEDB)
                messages.addWarningMessage(SQLITEDB)

            # otherwise put the database in the workspace
            else:
                messages.addMessage("Creating sqlite database in workspace")
                SQLITEDB = os.path.join(WORKSPACE, SQLITEDB)
        else:
            messages.addMessage("Sqlite db full path provided by user")

        # == TIMESTAMP HANDLING == #
        if not STARTTIMESTAMP or not ENDTIMESTAMP:
            messages.addMessage("No timestamp filter active")
            STARTTIMESTAMP = None
            ENDTIMESTAMP = None

        # == OUTPUT FILE NAME/TYPE HANDLING == #
        # if the full filepath of the output isn't specified, put in workspace
        if not os.path.dirname(OUTPUTGRID):
            #if geodatabase, don't use .shp extension
            if os.path.splitext(WORKSPACE)[1].lower() == ".gdb":
                if os.path.splitext(OUTPUTGRID)[1] == '.shp':
                    OUTPUTGRID = OUTPUTGRID[:-4]
            OUTPUTGRID = os.path.join(WORKSPACE,OUTPUTGRID)

        if not os.path.dirname(OUTPUTFOOTPRINTS):
            #if geodatabase, don't use .shp extension
            if os.path.splitext(WORKSPACE)[1].lower() == ".gdb":
                if os.path.splitext(OUTPUTFOOTPRINTS)[1] == '.shp':
                    OUTPUTFOOTPRINTS = OUTPUTFOOTPRINTS[:-4]
            OUTPUTFOOTPRINTS = os.path.join(WORKSPACE,OUTPUTFOOTPRINTS)


        # if the output shapefiles have no path, put them in the workspace
        arcpy.env.overwriteOutput = 1
        arcpy.env.workspace = WORKSPACE

        for RASTER in RASTER_LIST:
            messages.addMessage("Processing: {}".format(RASTER))
            proc = RasterProcessor(SQLITEDB,
                                   USNGGRID,
                                   OUTPUTGRID,
                                   OUTPUTFOOTPRINTS,
                                   ZONEFIELD,
                                   RASTER,
                                   STARTTIMESTAMP,
                                   ENDTIMESTAMP,
                                   messages)
            proc.process_raster()

            # if last raster, create the output files
            if RASTER == RASTER_LIST[-1]:
                proc.create_output()

        return
