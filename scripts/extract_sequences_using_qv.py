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
from Bio.Seq import reverse_complement, translate as bio_translate

def generate_placeholders(l):
    placeholder= ['%s'] # use ? For SQLite. See DBAPI paramstyle.
    return ', '.join(placeholder * l)


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def extract_feature_sql(server, seqfeature_ids, type=['CDS', 'rRNA', 'tRNA'], qualifier=['ID','locus_tag'], translate=False):
    """raw sql extraction of fasta seqfeatures
    """
    for chunk in chunks(seqfeature_ids, 900):
        sql = "SELECT f.seqfeature_id AS gid, \
                    fl.strand,\
                    substring(s.seq, fl.start_pos, (fl.end_pos - fl.start_pos)+1) AS subseq\
               FROM seqfeature f,\
                    location fl,\
                    biosequence s\
               WHERE f.seqfeature_id IN ({}) AND\
                     fl.seqfeature_id = f.seqfeature_id AND\
                     f.bioentry_id = s.bioentry_id".format(generate_placeholders(len(chunk)))

        features = server.adaptor.execute_and_fetchall(sql, tuple(chunk) )
        results = {}
        for sfid, strand, seq in features:
            results[sfid] = (strand, seq)

        qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv, term t WHERE seqfeature_id IN ({}) AND t.term_id = qv.term_id'.format(generate_placeholders(len(chunk)))
        qv = {}
        for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, tuple(chunk)):
            try:
                qv[seqfeature_id][name] = value
            except KeyError:
                qv[seqfeature_id] = {}
                qv[seqfeature_id][name] = value

        for seqfeature_id, (strand, seq) in results.items():
            name = str(seqfeature_id)
            for q in qualifier:
                try:
                    name += ' ' + qv[seqfeature_id][q]

                except KeyError:
                    pass
            if strand == -1:
                seq = reverse_complement(results[seqfeature_id][1])

            if translate:
                seq = bio_translate(seq)
            try:
                name += ' ' + qv[seqfeature_id]['product']
            except KeyError:
                pass

            print(">{}\n{}".format(name, seq))


def get_seqfeature_by_qv(server, qualifier, value, biodb=None):
    ''' find all seqfeatures that have the given value for the qualifier
        returns a list of seqfeature_id
    '''
    sql = "SELECT qv.seqfeature_id FROM seqfeature_qualifier_value qv join seqfeature s using(seqfeature_id) join bioentry b using(bioentry_id) join biodatabase bd using(biodatabase_id) WHERE term_id = (SELECT term_id FROM term WHERE name = %s AND ontology_id = (SELECT ontology_id from ontology WHERE name = 'Annotation Tags')) AND value = %s"
    if biodb:
        sql += " AND bd.name = %s"
        return server.adaptor.execute_and_fetchall(sql, (qualifier, value, biodb))
    else:
        return server.adaptor.execute_and_fetchall(sql, (qualifier, value))

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    seqfeature_ids = get_seqfeature_by_qv(server, args.qualifier, args.value, args.database_name)

    if args.output_format == 'feat-prot':
        extract_feature_sql(server, seqfeature_ids, type=['CDS'], translate=True )
    elif args.output_format == 'feat-nucl':
        extract_feature_sql(server, seqfeature_ids )
    #for dbname in server:
    #    db = server[dbname]
    #    for dbid, taxid in dbids.items():
    #        try:
    #            dbrec = db[dbid]
    #            if 'feat' in args.output_format:
    #                extract_feature(dbrec, args.output_format, sys.stdout)
    #            else:
    #                SeqIO.write(dbrec, sys.stdout, args.output_format)
    #        except KeyError:
    #            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host',
            default=5432)
    parser.add_argument('-u', '--user', help='database user name',
            default=getuser())
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to',
            default='localhost')
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['fasta', 'gb', 'gff', 'feat-prot', 'feat-nucl'], default='fasta')
    parser.add_argument('qualifier', help='name of the qualifier', default=None)
    parser.add_argument('value', help='value to match on' )
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

