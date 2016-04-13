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

def extract_feature_sql(server, dbids, type=['CDS', 'rRNA', 'tRNA'], qualifier=['ID','locus_tag'], translate=False):
    """raw sql extraction of fasta seqfeatures
    """
    for chunk in chunks(dbids, 900):


        seqs = dict(server.adaptor.execute_and_fetchall(\
                'SELECT bioentry_id, seq FROM biosequence WHERE bioentry_id IN ({})'.format(\
                generate_placeholders(len(chunk))),\
                tuple(chunk)))

        feat_select_simple_sql = 'SELECT seqfeature_id, bioentry_id FROM seqfeature WHERE bioentry_id IN ({}) '.format(generate_placeholders(len(chunk)))
        term_qualifier_sql = 'AND type_term_id in (SELECT term_id FROM term WHERE name in ({}))'.format(generate_placeholders(len(type)))
        features =  dict(server.adaptor.execute_and_fetchall(feat_select_simple_sql + term_qualifier_sql, tuple(chunk + type)))
        for feat_chunk in chunks(features.keys(), 900):

            location_select_sql = 'SELECT seqfeature_id, strand, start_pos, end_pos FROM location WHERE seqfeature_id IN ({})'.format(generate_placeholders(len(feat_chunk)))
            qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv, term t WHERE seqfeature_id IN ({}) AND t.term_id = qv.term_id'.format(generate_placeholders(len(feat_chunk)))
            qv = {}
            for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, tuple(feat_chunk)):
                try:
                    qv[seqfeature_id][name] = value
                except KeyError:
                    qv[seqfeature_id] = {}
                    qv[seqfeature_id][name] = value

            for seqfeature_id, strand, start_pos, end_pos in server.adaptor.execute_and_fetchall(location_select_sql, tuple(feat_chunk)):
                name = str(seqfeature_id)
                for q in qualifier:
                    try:
                        name += ' ' + qv[seqfeature_id][q]

                    except KeyError:
                        pass
                seq = seqs[features[seqfeature_id]][start_pos-1:end_pos]
                if strand == -1:
                    seq = reverse_complement(seq)

                if translate:
                    seq = bio_translate(seq)
                try:
                    name += ' ' + qv[seqfeature_id]['product']
                except KeyError:
                    pass

                print(">{}\n{}".format(name, seq))



def extract_feature(dbrec, output_format, fp, wanted_types=['CDS','rRNA', 'tRNA'], id_tag=None):

    for feature in dbrec.features:
        # only print the proteins
        if feature.type not in wanted_types:
            continue

        feat_extract = feature.extract(dbrec.seq.toseq())
        if output_format == 'feat-prot':
            feat_extract = feat_extract.translate()

        if id_tag:
            try:
                seqid = feature.qualifiers[id_tag][0]
            except KeyError:
                print("WARNING: cannot find {} tag for seqfeature {}".format(id_tag, feature._seqfeature_id), file=sys.stderr)
        else:
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

        feat_extract = SeqRecord(feat_extract, id=str(seqid), description=description)
        SeqIO.write(feat_extract, fp, 'fasta')

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    tax_name = False
    try:
        ncbi_tax = int(args.taxid)
    except ValueError:
        tax_name = True

    if not tax_name:
        taxon_id_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.ncbi_taxon_id  = %s AND include.right_value = include.left_value + 1)"

        rows = server.adaptor.execute_and_fetchall(taxon_id_lookup_sql, (ncbi_tax,))
    else:
        taxon_name_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name "\
                "WHERE name like %s) AND include.right_value = include.left_value + 1)"
        rows = server.adaptor.execute_and_fetchall(taxon_name_lookup_sql, (args.taxid,))

    dbids = {}
    for row in rows:
        dbids[(row[0], row[2])] = row[1]
    files = {}
    taxid_to_dbids = {}
    if args.split_species:
        taxon_file_mapping = {}
        for k, v in dbids.items():
            tname = server.adaptor.execute_and_fetch_col0("SELECT name from taxon_name where taxon_id = %s and name_class = %s", (v,'scientific name'))[0]
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
        for taxid, dbid_list in taxid_to_dbids.items():
            for dbid, dbname in dbid_list:
                db = server[dbname]
                seq_rec = db[dbid]
                outdata.setdefault(taxid, []).append(seq_rec)

        for taxid, dbrecs in outdata.items():
            with open(files[taxid], 'w') as fp:
                if 'feat' in args.output_format:
                    for dbrec in dbrecs:
                        extract_feature(dbrec, args.output_format, fp)
                else:
                    SeqIO.write(dbid_s, fp, args.output_format)

    else:
        if args.output_format == 'feat-prot':
            extract_feature_sql(server, [x[0] for x in dbids.keys()],type=['CDS'], translate=True )
        elif args.output_format == 'feat-nucl':
            extract_feature_sql(server, [x[0] for x in dbids.keys()])
        else:
            for (dbid, dbname), taxid in dbids.items():
                db = server[dbname]
                try:
                    dbrec = db[dbid]
                    SeqIO.write(dbrec, sys.stdout, args.output_format)
                except KeyError:
                    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences that are from either an NCBI taxonomy ID or the complete name of a taxonomic rank",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('database', help='name of biosql database')
    #parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host',
            default=5432)
    parser.add_argument('-u', '--user', help='database user name',
            default=getuser())
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to',
            default='localhost')
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences', choices=['fasta', 'gb', 'feat-prot', 'feat-nucl'], default='fasta')
    parser.add_argument('taxid', help='supply a ncbi taxonomy id that will be extracted. If an integer is supplied it will be interpreted as an NCBI taxonomy id; otherwise it will be interpreted as part of a taxonomy name (e.g. Proteobacteria)', default=None)
    parser.add_argument('-s', '--split_species', help='when there are multiple species to be returned, split them into separate files, based on their name, instead of printing to stdout', default=False, action='store_true')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

