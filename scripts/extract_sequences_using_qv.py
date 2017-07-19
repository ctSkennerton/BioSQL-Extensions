#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import sqlite3
from getpass import getpass, getuser
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from common import generate_placeholders, chunks, standard_options, \
        get_seqfeature_ids_from_qv, extract_feature_sql



def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    seqfeature_ids = get_seqfeature_ids_from_qv(server, args.qualifier, args.value, args.database_name, fuzzy=args.fuzzy)

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
    parser.add_argument('qualifier', help='name of the qualifier', default=None)
    parser.add_argument('value', help='value to match on' )
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

