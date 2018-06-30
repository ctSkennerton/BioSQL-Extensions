# -*- coding: utf-8 -*-
"""Unit test package for biosqlx."""

import os
import tempfile
from BioSQL import BioSeqDatabase

def temp_db_filename():
    """Generate a temporary filename for SQLite database."""
    # In memory SQLite does not work with current test structure since the tests
    # expect databases to be retained between individual tests.
    # TESTDB = ':memory:'
    # Instead, we use (if we can) /dev/shm
    try:
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db", dir='/dev/shm')
    except OSError:
        # We can't use /dev/shm
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db")
    os.close(h)
    return test_db_fname

def connection_parameters(create=False):
    """Get info for connecting to a database.

    :param create: Create a new, empty database.
    :returns: tuple of database name, database driver, database user, database password, database host
    """
    if create:
        dbname = temp_db_filename()
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(driver='sqlite3', db=dbname)
        try:
            server.load_database_sql(os.path.join(os.path.dirname(__file__), 'biosqldb-sqlite.sql'))
            server.commit()
            server.close()
        except Exception:
            # Failed, but must close the handle...
            server.close()
            raise
    else:
        dbname = os.path.join(os.path.dirname(__file__), 'test.db')

    return dbname, 'sqlite3', '', '', ''
