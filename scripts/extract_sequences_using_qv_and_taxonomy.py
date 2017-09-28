#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
from getpass import getpass, getuser
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from common import generate_placeholders, chunks, standard_options, \
        get_seqfeature_ids_from_qv, extract_feature_sql, \
        get_bioentries_from_taxonomy, get_bioseqid_for_seqfeature

def filter_seqfeature_ids_by_taxonomy(server, seqfeature_ids, taxid):
    wanted_bioentries = get_bioentries_from_taxonomy(server, taxid)
    seqfeature_bioentries = get_bioseqid_for_seqfeature(server, seqfeature_ids)
    seqfeature_bioentry_map = {}
    for row in seqfeature_bioentries:
        try:
            seqfeature_bioentry_map[row[1]].append(row[2])
        except KeyError:
            seqfeature_bioentry_map[row[1]] = [row[2]]

    final_seqfeatures = []
    for row in wanted_bioentries:
        try:
            final_seqfeatures.extend(seqfeature_bioentry_map[row[0]])
        except KeyError:
            pass

    return final_seqfeatures


def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    initial_seqfeature_ids = get_seqfeature_ids_from_qv(server, args.qualifier, args.value, args.database_name, fuzzy=args.fuzzy)
    seqfeature_ids = filter_seqfeature_ids_by_taxonomy(server, initial_seqfeature_ids, args.taxid)

    if args.feature_type is not None:
        types = args.feature_type
    elif args.output_format == 'feat-prot':
        types = ['CDS']
    elif args.output_format == 'feat-nucl':
        types = ['CDS', 'rRNA', 'tRNA']

    if args.output_format == 'feat-prot':
        extract_feature_sql(server, seqfeature_ids, type=types, translate=True )
    elif args.output_format == 'feat-nucl':
        extract_feature_sql(server, seqfeature_ids, type=types)


if __name__ == "__main__":
    parser = standard_options()
    parser.add_argument('-D', '--database-name', help='limit the extracted sequences from this namespace', dest='database_name')
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['feat-prot', 'feat-nucl'], default='feat-prot')
    parser.add_argument('-t', '--feature-type', help='restrict the results to feature type e.g. rRNA, tRNA, CDS. This option can be specified multiple times for multiple types', default=None, action='append')
    parser.add_argument('-f', '--fuzzy', help='the value can be a partial match', default=False, action='store_true')
    parser.add_argument('qualifier', help='name of the qualifier')
    parser.add_argument('value', help='value to match on' )
    parser.add_argument('taxid', help='supply a ncbi taxonomy id that will be extracted. If an integer is supplied it will be interpreted as an NCBI taxonomy id; otherwise it will be interpreted as part of a taxonomy name (e.g. Proteobacteria)')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)


