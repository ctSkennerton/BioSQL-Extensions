#!/usr/bin/env python
from __future__ import print_function
import sys
from BioSQL import BioSeqDatabase
from BioSQL import Loader
from common import standard_options, get_seqfeature_id_from_qv
import csv

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


def parse_input(infile, key):
    mapping = {}
    with open(infile) as fp:
        reader = csv.DictReader(fp, delimiter="\t")
        for row in reader:
            mapping[(key, row[key])] = {}
            for k, v in row.items():
                if k != key:
                    try:
                        mapping[(key, row[key])][k].append(v)
                    except KeyError:
                        mapping[(key, row[key])][k] = [v]
    return mapping

def parse_gff(infile):
    from BCBio import GFF
    mapping = {}
    with open(infile) as fp:
        for rec in GFF.parse(fp):
            for feature in rec.features:
                try:
                    protein_id = feature.qualifiers['ID'][0]
                    mapping[('ID', protein_id)] = feature.qualifiers
                except KeyError:
                    try:
                        protein_id = feature.qualifiers['locus_tag'][0]
                        mapping[('locus_tag', protein_id)] = feature.qualifiers
                    except KeyError:
                        pass
    return mapping


def add_annotation(db, mapping, isSeqfeatureAlready=False, replace=False):
    # this may be a little tricky depending on how the database is set up
    # since a bioentry is equivelent to a genbank file but genbank files could
    # be created from a whole chromosome or from an individual protein.
    # If it is from a single protein then the protein ID will be the bioentry_id
    # but if it is from a whole genome then it will be a seqfeature_id
    db_loader = CustomDBLoader(db.adaptor, db.dbid, False)

    for (term_name, protein), values in mapping.items():
        # Start by looking for bioentries that have the name
        if not isSeqfeatureAlready:
            seqfeature_id = get_seqfeature_id_from_qv(db, term_name, protein)

        else:
            seqfeature_id = int(protein)
        #print("loading ",kegg_id," into ",protein,"(",seqfeature_id,")")
        # now add in our qualifier and value onto that seqfeature
        if replace:
            for qualifier in values.keys():
                if qualifier == 'db_xref':
                    print('Cannot remove any current db_xref, you must do this manually for seqfeature {}'.format(seqfeature_id), file=sys.stderr)
                else:
                    #print("removing current annotations for <<{}>> tag in seqfeature {}".format(qualifier, seqfeature_id), file=sys.stderr)
                    db.adaptor.execute("delete from seqfeature_qualifier_value where term_id = \
                                (select term_id from term where ontology_id = \
                                    (select ontology_id from ontology where name = 'Annotation Tags')\
                                and name = %s)\
                            and seqfeature_id = %s", (qualifier, seqfeature_id))

        db_loader._load_seqfeature_qualifiers(values, seqfeature_id)

if __name__ == '__main__':
    from getpass import getpass
    parser = standard_options()
    parser.add_argument('-D', '--database-name', default=None,
        dest='dbname', help='name of the sub-database')

    # only one of the input formats can be given, but at least one must be given
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input', help='provide text file, tab delimited, where the first column is the name of the sequence feature and the following columns are the values of the annotation that you want to add. The first line must be a header line, which will be used as the name of the qualifier for the seqfeature.')
    group.add_argument('-g', '--gff', help='provide a gff3 formatted file whose attributes will be added to existing sequence features. Only the information in the last column of the gff file will be utilized so you must make sure that either the ID or locus_tag qualifiers are present in the gff file. If both are present then ID will be preferred over locus_tag. If neither are present then the record will be skipped. Make sure that the ID or locus_tag are unique (and present) in the database otherwise the attributes will not be loaded.')

    parser.add_argument('-s', '--seqfeature', help='The first column of the input file is the seqfeature id used by the database. Does not apply when using a gff file as input', action='store_true', default=False)
    parser.add_argument('--replace', help='replace any existing annotations for the given qualifiers', action='store_true', default=False)
    parser.add_argument('--key', help='name of the column that contains a unique identifier for the seqfeature. e.g. locus_tag', required=True)
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)

    server = BioSeqDatabase.open_database(driver=args.driver,
            db=args.database,
            user=args.user,
            host=args.host,
            passwd=args.password)

    if args.dbname is None:
        db = server[list(server.keys())[0]]
    else:
        db = server[args.dbname]

    if args.input is not None:
        mapping = parse_input(args.input, args.key)
    else:
        mapping = parse_gff(args.gff)

    add_annotation(db, mapping, args.seqfeature, args.replace)
    server.commit()

