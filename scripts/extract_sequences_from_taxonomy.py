#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import sqlite3
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_feature(dbrec, output_format, fp):

    for feature in dbrec.features:
        feat_extract = feature.extract(dbrec.seq.toseq())
        if output_format == 'feat-prot':
            feat_extract = feat_extract.translate()
        try:
            seqid = feature.qualifiers['ID'][0]
        except KeyError:
            try:
                seqid = feature.qualifiers['locus_tag'][0]
            except KeyError:
                seqid = feature._seqfeature_id

        description = ''
        try:
            description += ' '.join(feature.qualifiers['product'])
        except KeyError:
            pass
        else:
            description += ' '

        try:
            description += ' '.join(feature.qualifiers['gene'])
        except KeyError:
            pass

        feat_extract = SeqRecord(feat_extract, id=seqid, description=description)
        SeqIO.write(feat_extract, fp, 'fasta')

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    try:
        ncbi_tax = int(args.taxid)
        taxon_id_lookup_sql = "SELECT bioentry_id, taxon_id FROM bioentry WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.ncbi_taxon_id  = %s) AND include.right_value = include.left_value + 1)"

        rows = server.adaptor.execute_and_fetchall(taxon_id_lookup_sql, (ncbi_tax,))
    except:
        taxon_name_lookup_sql = "SELECT bioentry_id, taxon_id FROM bioentry WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name "\
                "WHERE name like %s) AND include.right_value = include.left_value + 1)"

        rows = server.adaptor.execute_and_fetchall(taxon_name_lookup_sql, (args.taxid,))
    finally:
        dbids = dict(rows)
        files = {}
        taxid_to_dbids = {}
        if args.split_species:
            taxon_file_mapping = {}
            for k, v in dbids.items():
                tname = server.adaptor.execute_and_fetch_col0("SELECT name from taxon_name where taxon_id = %s", (v,))[0]
                tname = tname.replace(' ', '_')
                if args.output_format == 'gb':
                    tname += '.gb'
                elif args.output_format == 'feat-prot':
                    tname += '.faa'
                else:
                    tname += '.fna'
                files[v] = tname
                taxid_to_dbids.setdefault(v, []).append(k)


    if args.split_species:
        # got to save all of the records before printing them out
        outdata = {}
        for db in server.values():
            for taxid, dbid_list in taxid_to_dbids.items():
                for dbid in dbid_list:
                    try:
                        seq_rec = db[dbid]
                        outdata.setdefault(taxid, []).append(seq_rec)
                    except KeyError:
                        pass

        for taxid, dbrecs in outdata.items():
            with open(files[taxid], 'w') as fp:
                if 'feat' in args.output_format:
                    for dbrec in dbrecs:
                        extract_feature(dbrec, args.output_format, fp)
                else:
                    SeqIO.write(dbid_s, fp, args.output_format)

    else:
        for dbname in server:
            db = server[dbname]
            for dbid, taxid in dbids.items():
                try:
                    dbrec = db[dbid]
                    if 'feat' in args.output_format:
                        extract_feature(dbrec, args.output_format, sys.stdout)
                    else:
                        SeqIO.write(dbrec, sys.stdout, args.output_format)
                except KeyError:
                    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host')
    parser.add_argument('-u', '--user', help='database user name')
    parser.add_argument('-P', '--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to', default='localhost')
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['fasta', 'gb', 'feat-prot', 'feat-nucl'], default='fasta')
    parser.add_argument('taxid', help='supply a ncbi taxonomy id that will be extracted. If an integer is supplied it will be interpreted as an NCBI taxonomy id; otherwise it will be interpreted as part of a taxonomy name (e.g. Proteobacteria)', default=None)
    parser.add_argument('-s', '--split_species', help='when there are multiple species to be returned, split them into separate files, based on their name, instead of printing to stdout', default=False, action='store_true')
    args = parser.parse_args()
    main(args)

