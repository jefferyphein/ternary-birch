import sqlite3
import operator
import array
import logging

from functools import reduce

from sage.all import vector as sage_vector
from sage.all import Integers
from sage.all import Integer

class EigenvectorDatabase:
    ###########################################################################
    # SQL instructions for building a database from scratch.

    CREATE_TABLE_GENUS = '''CREATE TABLE IF NOT EXISTS genus (
        level               INTEGER     PRIMARY KEY,
        ramified_primes     TEXT        NOT NULL,
        num_eigenvectors    INTEGER     NOT NULL DEFAULT -1,
        seed                INTEGER     NOT NULL
    )'''

    CREATE_TABLE_EIGENVECTORS = '''CREATE TABLE IF NOT EXISTS eigenvectors (
        id                  INTEGER     PRIMARY KEY,
        level               INTEGER     NOT NULL,
        conductor           INTEGER     NOT NULL,
        vector              TEXT        NOT NULL,
        FOREIGN KEY (level) REFERENCES genus (level)
    )'''

    CREATE_TABLE_EIGENVALUES = '''CREATE TABLE IF NOT EXISTS eigenvalues (
        vector_id           INTEGER     NOT NULL,
        p                   INTEGER     NOT NULL,
        value               INTEGER     NOT NULL,
        UNIQUE (vector_id, p),
        FOREIGN KEY (vector_id) REFERENCES eigenvectors (id)
    )'''

    CREATE_INDEX_EIGENVECTORS = '''CREATE INDEX IF NOT EXISTS eigenvalues_index
        ON eigenvectors (level, conductor)'''

    CREATE_INDEX_EIGENVALUES = '''CREATE INDEX IF NOT EXISTS eigenvectors_index
        ON eigenvalues (vector_id, p)'''

    ###########################################################################
    # Subroutines required to make this class a context manager.

    def __init__(self, filename):
        self.filename = filename
        self.conn = self.create_connection(self.filename)
        self.conn.row_factory = sqlite3.Row
        self.create_database(self.conn)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.conn.close()

    ###########################################################################
    # Build the genus objects in a variety of ways.

    def get_genus(self, level):
        cur = self.conn.cursor()

        cur.execute("SELECT * FROM genus WHERE level=?", (int(level),))
        result = cur.fetchone()

        if result is None:
            # If no results are found, construct a new genus.
            genus = BirchGenus(level)

            ramified_primes = EigenvectorDatabase.pack_ramified_primes(genus)
            values = (int(level), ramified_primes, int(genus.seed()))
            cur.execute("INSERT INTO genus (level, ramified_primes, seed) VALUES (?, ?, ?)", values)
        else:
            # Otherwise, construct the genus from stored parameters.
            assert level == result['level']
            ramified_primes = EigenvectorDatabase.unpack_ramified_primes(result['ramified_primes'])
            seed = result['seed']
            genus = BirchGenus(level, ramified_primes=ramified_primes, seed=seed)

        cur.close()

        return genus

    def get_genus_with_eigenvectors(self, level, precise=True, remove_oldforms=True):
        genus = self.get_genus(level)
        level = int(level)

        cur = self.conn.cursor()

        cur.execute("SELECT num_eigenvectors FROM genus WHERE level=?", (level,))
        result = cur.fetchone()
        num_eigenvectors = result['num_eigenvectors']
        if num_eigenvectors < 0:
            # We haven't yet computed the rational eigenvectors, so let's do
            # that now.
            evecs = genus.rational_eigenvectors(precise=precise)

            # Compute eigenvalues up to a specific bound so that we can
            # distinguish oldforms from newforms.
            # TODO: This should probably be the Sturm bound???
            # NOTE: Probably doesn't matter, as long as we compute *enough*
            #   eigenvalues to allow us to easily identify the duplicates that
            #   we know will be present.
            genus.compute_eigenvalues_upto(100, precise=precise)

            # If there are an even number of ramified primes, we need to remove
            # oldforms.
            if remove_oldforms:
                # Due to the way that we have chosen ramified primes, oldforms
                # will only occur when the level (assumed squarefree) is the
                # product of an even number of primes. In this case, we have
                # assumed that all but the largest prime factor is ramified,
                # and so oldforms will occur exactly twice to account for the
                # additional Atkin-Lehner sign we've introduced at the
                # unramified prime.
                #
                # For example, at level p*q with p<q, oldforms arise from
                # level p with the A-L sign inherited from level p. A newform
                # at level p with - sign will then occur at level p*q with
                # signs -+ and -- at level p*q.
                #
                # This means that we can easily identify oldforms due to their
                # multiplicity and remove them directly without having to check
                # for duplicates at level p.

                # Determine the oldform level.
                old_level = reduce(operator.mul, genus.ramified_primes())

                # If the old level and the current new level match, do nothing.
                if old_level != level:
                    # Determine the large unramified prime.
                    q = level // old_level

                    # Identify eigenvectors with duplicate eigenvalues.
                    removal = sorted(EigenvectorDatabase.identify_duplicates(evecs, q))

                    # Remove oldforms, starting at the end of the list.
                    for n in removal[::-1]:
                        del evecs[n]

                    # Since we may have removed some eigenvectors, we must
                    # reset the eigenvector manager before we can compute
                    # eigenvalues again.
                    if len(removal) > 0:
                        genus.reset_eigenvector_manager()

            try:
                # Loop over each eigenvector and create an entry in the database
                # for each one.
                for vec in evecs:
                    coords = EigenvectorDatabase.pack_vector(vec['vector'])
                    values = (level, int(vec['conductor']), coords)
                    cur.execute("INSERT INTO eigenvectors (level, conductor, vector) VALUES (?, ?, ?)", values)

                    # Recover the newly assigned vector id so we can insert its
                    # known eigenvalues.
                    vector_id = cur.lastrowid
                    vec['vector_id'] = vector_id

                    # Build a list of (vector_id, prime, eigenvalue) tuples to be
                    # inserted into the eigenvalue table.
                    eigenvalues = []
                    for item in vec['aps'].items():
                        eigenvalues.append((vector_id, int(item[0]), int(item[1])))

                    cur.executemany("INSERT INTO eigenvalues (vector_id, p, value) VALUES (?, ?, ?)", eigenvalues)

                # Now that all of the eigenvectors are in the table, update the
                # corresponding entry to both indicate that we've computed
                # eigenvectors, but also how many we obtained.
                cur.execute("UPDATE genus SET num_eigenvectors=? WHERE level=?", (len(evecs), level))

            except KeyboardInterrupt as e:
                self.conn.rollback()
                raise KeyboardInterrupt("Rolling back database to previous commit.")

        elif num_eigenvectors > 0:
            # Select all of the eigenvectors from the database.
            cur.execute("SELECT * FROM eigenvectors WHERE level=?", (level,))
            results = cur.fetchall()
            assert len(results) > 0

            # Construct the vector objects to be added to the genus.
            vector_ids = []
            vectors = dict()
            for result in results:
                vector_id = result['id']
                vector_ids.append(vector_id)
                vec = dict()
                vec['conductor'] = Integer(result['conductor'])
                vec['vector'] = EigenvectorDatabase.unpack_vector(result['vector'])
                vec['aps'] = dict()
                vec['vector_id'] = vector_id
                vectors[vector_id] = vec

            # Select all eigenvalues from the database for each eigenvector.
            cur.execute("SELECT * FROM eigenvalues WHERE vector_id IN ({0})".format(", ".join("?" for id in vector_ids)), vector_ids)
            results = cur.fetchall()
            for result in results:
                vector_id = result['vector_id']
                p = Integer(result['p'])
                value = Integer(result['value'])
                vectors[vector_id]['aps'][p] = value

            # Add each eigenvector to the genus object.
            for vec in vectors.items():
                genus.add_eigenvector(vec[1])

        cur.close()

        return genus

    def update_database_from_genus(self, genus, precise=True):
        level = int(genus.level())
        seed = int(genus.seed())
        ramified_primes = EigenvectorDatabase.pack_ramified_primes(genus)

        cur = self.conn.cursor()

        # Query the database to make sure this genus is present.
        values = (level, seed, ramified_primes)
        cur.execute("SELECT * FROM genus WHERE level=? AND seed=? AND ramified_primes=?", values)
        result = cur.fetchone()

        if result is None:
            raise Exception("No genus found in the database associated to this level, seed, and set of ramified primes.")

        # Do not proceed if we don't know how many eigenvectors exist.
        if result['num_eigenvectors'] < 0:
            raise Exception("Unknown number of eigenvectors associated to this genus.")

        # Nothing to do, so return.
        if result['num_eigenvectors'] == 0:
            return

        # Recover the eigenvectors from the genus. If used correctly, this will
        # return the already existing list of eigenvectors.
        eigenvectors = genus.rational_eigenvectors(precise=precise)

        try:
            for vec in eigenvectors:
                if not 'vector_id' in vec:
                    raise Exception("No vector id associated to eigenvector within genus... TODO: Resolve this automatically instead of throwing an exception.")
                vector_id = vec['vector_id']

                # Accumulate all (vector_id, prime, eigenvalue) tuples into a list
                # and then add each of them to the database. Conflicts are silently
                # ignored by the database.
                eigenvalues = []
                for p,value in vec['aps'].items():
                    values = (vector_id, int(p), int(value))
                    eigenvalues.append(values)
                cur.executemany("INSERT OR IGNORE INTO eigenvalues (vector_id, p, value) VALUES (?, ?, ?)", eigenvalues)
        except KeyboardInterrupt as e:
            self.conn.rollback()
            raise KeyboardInterrupt("Rolling back database to previous commit.")

        cur.close()

    ###########################################################################
    # Useful subroutines for managing when commits are flushed into the
    # database as well as gracefully shutting down and rolling back the
    # database.

    def commit(self):
        self.conn.commit()

    def rollback(self):
        self.conn.rollback()

    def close(self):
        self.conn.close()

    ###########################################################################
    # Static helper methods used for setting up and connecting to the database.

    @staticmethod
    def create_connection(filename):
        return sqlite3.connect(filename)

    @staticmethod
    def create_database(conn):
        c = conn.cursor()
        c.execute("PRAGMA foreign_keys = 1")
        c.execute(EigenvectorDatabase.CREATE_TABLE_GENUS)
        c.execute(EigenvectorDatabase.CREATE_TABLE_EIGENVECTORS)
        c.execute(EigenvectorDatabase.CREATE_TABLE_EIGENVALUES)
        c.execute(EigenvectorDatabase.CREATE_INDEX_EIGENVECTORS)
        c.execute(EigenvectorDatabase.CREATE_INDEX_EIGENVALUES)

    ###########################################################################
    # Static helper function for identifying eigenvectors with duplicate
    # eigenvectors.

    @staticmethod
    def identify_duplicates(vectors, prime):
        removal = []
        for n,vec1 in enumerate(vectors):
            cond = vec1['conductor']
            if cond % prime != 0: continue
            for m,vec2 in enumerate(vectors):
                if vec2['conductor'] * prime != cond: continue
                match = True
                for p,value in vec1['aps'].items():
                    if match and p in vec2['aps'] and value != vec2['aps'][p]:
                        match = False
                if match:
                    removal.append(n)
                    removal.append(m)
        return removal

    ###########################################################################
    # Static helper methods used for packing/unpacking data to/from database.

    @staticmethod
    def pack_ramified_primes(genus):
        return ",".join(map(str, genus.ramified_primes()))

    @staticmethod
    def unpack_ramified_primes(s):
        return list(map(Integer, s.split(",")))

    @staticmethod
    def pack_vector(vec):
        return ",".join(map(str, list(vec)))

    @staticmethod
    def unpack_vector(s):
        return sage_vector(Integers(), map(int, s.split(",")))

###############################################################################
# Do stuff.

def main():
    logger = logging.getLogger('compute_eigenvalues')
    logger.setLevel(logging.INFO)

    handler = logging.FileHandler('output.log')
    handler.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    with EigenvectorDatabase("eigenvectors.sqlite3") as db:
        last_commit = datetime.now()
    
        for level in range(2, 10000):
            # We can only recover newforms at squarefree levels.
            if not Integer(level).is_squarefree(): continue
    
            # Do the computations to determine the newforms, as well as their
            # Hecke eigenvalues for p < 1000.
            start_time = datetime.now()
            genus = db.get_genus_with_eigenvectors(level, remove_oldforms=True, precise=False)
            eigenvectors = genus.rational_eigenvectors()
            genus.compute_eigenvalues_upto(1000, precise=False)
            db.update_database_from_genus(genus, precise=False)
            end_time = datetime.now()
    
            logger.info("level = %-7s  dim = %-6s  num_vecs = %-3s  time = %s",
                level, genus.dimension(), len(eigenvectors), end_time-start_time)
    
            # Commit changes to the database every minute, so that if we cancel,
            # we don't lose too much data.
            if (end_time-last_commit).total_seconds() >= 60.0:
                print("Checkpoint at {}".format(level))
                db.commit()
                last_commit = end_time

        # Commit any final changes before we exit the context manager.
        db.commit()

if __name__ == "__main__":
    main()
