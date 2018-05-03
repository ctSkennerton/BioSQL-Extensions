#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
import csv
from getpass import getpass
from common import standard_options, generate_placeholders, chunks


def dbxref_dict(server, seqfeature_ids):
    db_qv = {}
    for feat_chunk in chunks(seqfeature_ids, 900):
        #sql = "SELECT s.seqfeature_id, d.dbname || ':' || d.accession AS kegg_id "\
        #        "FROM seqfeature_dbxref s "\
        #        "JOIN dbxref d USING(dbxref_id) "\
        #        "WHERE s.seqfeature_id IN ({})".format(generate_placeholders(len(feat_chunk)))
        sql = "SELECT s.seqfeature_id, d.dbname, d.accession, t.name, dqv.value "\
                "FROM seqfeature_dbxref s "\
                "JOIN dbxref d USING(dbxref_id) "\
                "LEFT JOIN dbxref_qualifier_value dqv USING(dbxref_id) "\
                "LEFT JOIN term t USING(term_id) "\
                "WHERE s.seqfeature_id IN ({}) "\
                "ORDER BY s.seqfeature_id, d.dbname, s.rank".format(generate_placeholders(len(feat_chunk)))
        for seqfeature_id, dbname, dbxref, name, value in server.adaptor.execute_and_fetchall(sql, tuple(feat_chunk)):
        #for seqfeature_id, dbxref in server.adaptor.execute_and_fetchall(sql, tuple(feat_chunk)):
            try:
                db_qv[seqfeature_id][dbname] = dbxref
            except KeyError:
                db_qv[seqfeature_id] = {}
                db_qv[seqfeature_id][dbname] = dbxref

            if name:
                db_qv[seqfeature_id][name] = value
    return db_qv



def qv_dict(server, seqfeature_ids):
    qv = {}
    for feat_chunk in chunks(seqfeature_ids, 900):
        feat_chunk2 = tuple(feat_chunk)
        qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv JOIN term t ON t.term_id = qv.term_id WHERE seqfeature_id IN ({})'.format(generate_placeholders(len(feat_chunk)))

        taxonomy_sql = 'SELECT seqfeature_id, bioentry.name, biodatabase.name, lineage.lineage FROM seqfeature JOIN bioentry USING(bioentry_id) JOIN biodatabase USING(biodatabase_id) LEFT JOIN lineage ON taxon_id = lineage.id WHERE seqfeature_id IN ({})'.format(generate_placeholders(len(feat_chunk)))
        for seqfeature_id, contig, namespace, lineage in server.adaptor.execute_and_fetchall(taxonomy_sql, feat_chunk2):
            try:
                if lineage:
                    qv[seqfeature_id]['taxonomy'] = lineage
                    qv[seqfeature_id]['organism'] = lineage.split(';')[-1]

                qv[seqfeature_id]['bioentry'] = contig
                qv[seqfeature_id]['sample'] = namespace
            except KeyError:
                qv[seqfeature_id] = {}
                if lineage:
                    qv[seqfeature_id]['taxonomy'] = lineage
                    qv[seqfeature_id]['organism'] = lineage.split(';')[-1]

                qv[seqfeature_id]['bioentry'] = contig
                qv[seqfeature_id]['sample'] = namespace

        for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, feat_chunk2):
            if not name:
                continue
            try:
                qv[seqfeature_id][name] = value
            except KeyError:
                qv[seqfeature_id] = {}
                qv[seqfeature_id][name] = value

    return qv

def print_feature_qv_csv(server, sfids, outfile=sys.stdout):
    """raw sql extraction of fasta seqfeatures
    """
    qv = qv_dict(server, sfids)
    dbxr = dbxref_dict(server, sfids)
    for sf, data in dbxr.items():
        for n, v in data.items():
            qv[sf][n] = v

    columns = set()
    for sf, data in qv.items():
        columns |= set(data.keys())

    columns = sorted(columns)
    writer = csv.writer(outfile)
    writer.writerow(['seqfeature'] + columns)
    for sf, data in qv.items():
        row = [sf]
        for i in columns:
            row.append(data.get(i, None))
        writer.writerow(row)


def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)
    sfids = []
    with open(args.input) as fp:
        for line in fp:
            sfids.append(line.rstrip())

    print_feature_qv_csv(server, sfids)


if __name__ == "__main__":
    parser = standard_options()
    parser.add_argument('input', help='file containing seqfeature ids')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

