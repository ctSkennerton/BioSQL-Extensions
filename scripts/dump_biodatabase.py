#!/usr/bin/env python
import sys
from getpass import getpass
from BioSQL import BioSeqDatabase
from common import standard_options, generate_placeholders, chunks, extract_feature_sql

def get_seqfeature_for_db(server, biodb):
    ''' find all seqfeatures that have the given value for the qualifier
        returns a list of seqfeature_id
    '''
    sql = "SELECT qv.seqfeature_id FROM seqfeature_qualifier_value qv join seqfeature s using(seqfeature_id) join bioentry b using(bioentry_id) join biodatabase bd using(biodatabase_id) WHERE bd.name = %s"
    return server.adaptor.execute_and_fetchall(sql, (biodb,))

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    if args.output_format == 'fasta':
        from Bio import SeqIO
        db = server[args.database_name]
        for rec in db.values():
            SeqIO.write(rec, sys.stdout, args.output_format)
    else:

        seqfeature_ids = get_seqfeature_for_db(server, args.database_name)

        if args.output_format == 'feat-prot':
            extract_feature_sql(server, seqfeature_ids, type=['CDS'], translate=True )
        elif args.output_format == 'feat-nucl':
            extract_feature_sql(server, seqfeature_ids )


if __name__ == "__main__":
    parser = standard_options()
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', required=True)
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['feat-prot', 'feat-nucl', 'fasta'], default='feat-prot')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

