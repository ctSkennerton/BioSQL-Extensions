#!/usr/bin/env python
from __future__ import print_function
import sys
import sqlite3
from BioSQL import BioSeqDatabase
from BioSQL import Loader

def parse_kegg(infile):
    mapping = {}
    for line in infile:
        try:
            protein, kegg = line.rstrip().split('\t')
            mapping[protein] = kegg
        except:
            pass

    return mapping

def add_kegg_annotation(db, mapping):
    # this may be a little tricky depending on how the database is set up
    # since a bioentry is equivelent to a genbank file but genbank files could
    # be created from a whole chromosome or from an individual protein.
    # If it is from a single protein then the protein ID will be the bioentry_id
    # but if it is from a whole genome then it will be a seqfeature_id

    db_loader = Loader.DatabaseLoader(db.adaptor, db.dbid, False)

    # We need the internal ID for a CDS and locus_tag types for later
    cds_term_id = db.adaptor.execute_and_fetchall('select term_id from term where name = "CDS"')[0][0]
    locus_tag_term_id = db.adaptor.execute_and_fetchall('select term_id from term where name = "locus_tag"')[0][0]
    for protein, kegg_id in mapping.items():
        # Start by looking for bioentries that have the name
        try:
            seqid = db.adaptor.fetch_seqid_by_display_id(db.dbid, protein)

            # ok now lets add in our new qualifier to to the CDS feature of this bioentry
            # start by finding the seqfeature id so that we can add qualifiers
            sql = r'select sf.seqfeature_id from seqfeature sf where bioentry_id = ? and sf.type_term_id = ?;'
            seqfeature_id = db.adaptor.execute_and_fetchall(sql, (seqid, cds_term_id))[0][0]

        except IndexError, e:
            #print(e, file=sys.stderr)

            # so if that fails then look for the name in the locus tag of the seqfeatures
            sql = r'select seqfeature_id from seqfeature_qualifier_value where term_id = ? and value = ?'
            try:
                seqfeature_id = db.adaptor.execute_and_fetchall(sql, (locus_tag_term_id, protein))[0][0]
            except IndexError, e:
                print(e, file=sys.stderr)
                raise e

        print("loading ",kegg_id," into ",protein,"(",seqfeature_id,")")
        # now add in our qualifier and value onto that seqfeature
        db_loader._load_seqfeature_qualifiers({'db_xref': ["KO:"+kegg_id]}, seqfeature_id)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', help='pre-created biosql database')
    parser.add_argument('-D', '--database-name', default='metagenomic_database',
        dest='dbname', help='name of the sub-database')
    parser.add_argument('-k', '--kegg', type=argparse.FileType(),
        help='input file containing a two column, tab delimited, file where '\
        'the first column is the nume of a protein stored in the database '\
        'and the second column is the KEGG ontology of the protein')
    args = parser.parse_args()

    mapping = parse_kegg(args.kegg)
    server = BioSeqDatabase.open_database(driver='sqlite3',db=args.database)
    db = server[args.dbname]
    add_kegg_annotation(db, mapping)
    server.commit()

