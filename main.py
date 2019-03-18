from __future__ import with_statement
import os
import arcpy
import glob
import uuid
import string
import random
import hashlib
import sqlite3
from arcpy import env
from arcpy.sa import SetNull

from contextlib import closing


# Globals
# TODO Set as tool/gui inputs
ROOTDIR = r'G:\ssdworkspace\raster_to_usng_rollup\data'
WORKSPACE = r'{}\scratch'.format(ROOTDIR)
SQLITEDB = r'{}\db\rollupdata.sqlite'.format(ROOTDIR)
USNGGRID = r'{}\testdata\usng_1k_withpr.shp'.format(ROOTDIR)
ZONEFIELD = "USNG_1KM"
RASTER = r'{}\testdata\Irma_DG_OB_FEMA.tif'.format(ROOTDIR)
#DISCARDCELLS = "Value <= 0"
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
            arcpy.gp.CreateSQLiteDB(self.db, 'SPATIALITE')

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
    def __init__(self, db, poly, zonefield, raster, discard_condition = None):
        self.db = db
        self.poly = poly
        self.polyname = os.path.basename(poly)
        self.polyhash = self.generate_md5(poly) #TODO generate unique from GDB
        self.zonefield = zonefield
        self.raster = raster
        self.rastername = os.path.basename(raster)
        self.rasterhash = self.generate_md5(raster) #TODO generate unique from GDB
        self.discard_condition = discard_condition
        self.null_raster = None
        self.zonal_stats_keys = None
        self.zonal_stats_data = None

    def process_raster(self):
        if self.precheck_db():
            print("Warning: Identical data found in database")
            print("Discontinuing processing")
            return

        if self.discard_condition:
            self.process_null_raster()
            self.get_zonal_stats_np_array(self.null_raster)
        else:
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
        #self.load_footprint_table(rasterid)

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
                  discard_condition TEXT);""")

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
                  the_geom POLYGON,
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
                INSERT into zonal_stats(raster_id, feature_id, count, min, max, mean, std)
                VALUES (?,?,?,?,?,?,?)
            """
            process_list = [[rasterid] + list(row) for row in self.zonal_stats_data]
            dbc.executemany(insert_sql, (process_list))

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
                SELECT r.discard_condition
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
                                  'all')

            # check if the discard conditions has been used before on these data
            if results:
                for result in results:
                    if result[0] == self.discard_condition:
                        return True
                return False
            else:
                return False

    def load_raster_table(self):
        with DBC(self.db) as dbc:
            insert_sql = """
                INSERT INTO raster(name, hash, discard_condition) 
                VALUES (?,?,?)
                """
            rowid = dbc.execute(insert_sql, (self.rastername, self.rasterhash, self.discard_condition), 'rowid')
            select_sql = """
                SELECT id FROM raster WHERE rowid = ?
                """
            rasterid = dbc.execute(select_sql, (rowid), 'one')

            return rasterid[0]

    def generate_md5(self, file_name):
        """
        Generate the md5sum of an input file without reading it all to memory

        :param file_name: Input file to be hashed, must be readable
        :type file_name: String
        :return: Hash value for file at location
        :rtype: String
        """
        hash_md5 = hashlib.md5()
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
        no_data_val = raster_desc.noDataValue
        null_raster = SetNull(self.raster, self.raster, self.discard_condition)

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

    def get_zonal_stats_np_array(self, raster):
        """
        Create the raster stats data by a field in the polygon

        :param raster: Raster location
        :type raster: String
        """

        # create a 10 digit random filename with a non numeric first
        random_letter = random.choice(string.ascii_lowercase)
        random_uuid = str(uuid.uuid4().hex)[:9]
        temp_file = "{}{}.dbf".format(random_letter, random_uuid)

        temp_dbf = os.path.join(arcpy.env.workspace, temp_file)
        # Run the zonal statistics to table
        arcpy.gp.ZonalStatisticsAsTable_sa(self.poly, self.zonefield,
                                           raster, temp_dbf,
                                           "DATA", "ALL")

        zonal_stats_data = arcpy.da.TableToNumPyArray(
            temp_dbf,[self.zonefield, "COUNT", "MIN", "MAX", "MEAN", "STD"])
        zonal_stats_keys = ["zone", "count", "min", "max", "mean", "std"]

        # remove intermediate data from filesystem
        self.arcpy_rm_file(temp_dbf+'*', glob_opt=True)

        self.zonal_stats_keys = zonal_stats_keys
        self.zonal_stats_data = zonal_stats_data

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


if __name__ == "__main__":

    # Config
    env.overwriteOutput = 1
    env.workspace = WORKSPACE

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    proc = RasterProcessor(SQLITEDB, USNGGRID, ZONEFIELD, RASTER, DISCARDCELLS)
    #proc.clear_db()
    proc.truncate_db()
    proc.process_raster()
