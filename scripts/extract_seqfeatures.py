#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
from getpass import getpass, getuser
from BioSQL import BioSeqDatabase
from common import standard_options, extract_feature_sql

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)
    seqfeature_ids = []
    with open(args.infile) as fp:
        for line in fp:
            seqfeature_ids.append(int(line.rstrip()))

    if args.output_format == 'feat-prot':
        extract_feature_sql(server, seqfeature_ids, type=['CDS'], translate=True )
    elif args.output_format == 'feat-nucl':
        extract_feature_sql(server, seqfeature_ids )


if __name__ == "__main__":
    from getpass import getpass
    parser = standard_options()
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['feat-prot', 'feat-nucl'], default='feat-prot')
    parser.add_argument('infile', help='file containing seqfeature_ids, one per line')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

