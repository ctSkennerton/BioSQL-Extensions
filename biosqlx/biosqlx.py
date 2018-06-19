# -*- coding: utf-8 -*-

from BioSQL import BioSeqDatabase
from BioSQL import Loader
class CustomDBLoader(Loader.DatabaseLoader):
    """This is a slightly modified version of Loader.DatabaseLoader
    """
    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
            """Insert the (key, value) pair qualifiers relating to a feature (PRIVATE).
            Qualifiers should be a dictionary of the form:
                {key : [value1, value2]}
                Before insertion, each qualifier will be checked to make sure that there
                is no collision between current annotations, if there is the rank of the
                new features will be changed to prevent collisions
            """
            tag_ontology_id = self._get_ontology_id('Annotation Tags')
            for qualifier_key in qualifiers:
                # Treat db_xref qualifiers differently to sequence annotation
                # qualifiers by populating the seqfeature_dbxref and dbxref
                # tables.  Other qualifiers go into the seqfeature_qualifier_value
                # and (if new) term tables.
                if qualifier_key != 'db_xref':
                    qualifier_key_id = self._get_term_id(qualifier_key,
                                                         ontology_id=tag_ontology_id)
                    # now add all of the values to their table
                    entries = qualifiers[qualifier_key]
                    if not isinstance(entries, list):
                        # Could be a plain string, or an int or a float.
                        # However, we exect a list of strings here.
                        entries = [entries]
                    # check for any rows for this term_id in this seqfeature_id
                    # and return the max rank
                    sql = "SELECT max(rank) FROM seqfeature_qualifier_value" \
                              " WHERE term_id = %s AND seqfeature_id = %s"
                    qual_rank_start = self.adaptor.execute_and_fetchall(sql, \
                            (qualifier_key_id, seqfeature_id))[0][0]
                    if not qual_rank_start:
                        # if there are no rows
                        qual_rank_start = 1

                    for i, qualifier_value in enumerate(entries):
                        # -1 cause enumerate begins with 1 rather than 0
                        qual_value_rank = i + qual_rank_start - 1
                        #print("qual_rank_start", qual_rank_start, i, qual_rank_start)
                        qualifier_value = entries[qual_value_rank]
                        if qualifier_value != "":
                            sql = r"INSERT INTO seqfeature_qualifier_value "\
                                  r" (seqfeature_id, term_id, rank, value) VALUES"\
                                  r" (%s, %s, %s, %s)"
                            self.adaptor.execute(sql, (seqfeature_id,
                                                       qualifier_key_id,
                                                       qual_value_rank,
                                                       qualifier_value))
                else:
                    # The dbxref_id qualifier/value sets go into the dbxref table
                    # as dbname, accession, version tuples, with dbxref.dbxref_id
                    # being automatically assigned, and into the seqfeature_dbxref
                    # table as seqfeature_id, dbxref_id, and rank tuples
                    self._load_seqfeature_dbxref(qualifiers[qualifier_key],
                                                 seqfeature_id)

