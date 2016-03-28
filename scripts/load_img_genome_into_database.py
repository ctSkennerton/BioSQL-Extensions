#!/usr/bin/env python
import sys
from os.path import split as path_split, join as path_join, isfile
import argparse
import sqlite3
import re
from tempfile import TemporaryFile
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio import Entrez
from BCBio import GFF
Entrez.email = 'c.skennerton@gmail.com'

def add_sequence_dbxref(inIter, taxid, img_genome_id):
    inIter.annotations['ncbi_taxid'] = taxid
    try:
        inIter.dbxrefs.append("IMG_genome:" + img_genome_id)
    except KeyError:
        inIter.dbxrefs = ["IMG_genome:" + img_genome_id]
    return inIter


def add_gene_dbxref(xref_file, kegg_file):

    ret = {}

    def add_or_append(key, value):
        try:
            ret[key].append(value)
        except KeyError:
            ret[key] = [value]


    def add_kegg(kegg_file):
        with open(kegg_file) as fp:
            # skip the first line since that is the header
            next(fp)
            for line in fp:
                fields = line.rstrip().split("\t")
                add_or_append(fields[0], 'ko:' + fields[9][3:])

    def add_xref(xref_file):
        with open(xref_file) as fp:
            # skip the first line since that is the header
            next(fp)
            for line in fp:
                fields = line.rstrip().split("\t")
                if 'GenBank' in fields[1]:
                    add_or_append(fields[0], 'GenBank:' + fields[2])
                elif 'UniProt' in fields[1]:
                    add_or_append(fields[0], 'UniProtKB:' + fields[2])

    if isfile(xref_file):
        add_xref(xref_file)

    if isfile(kegg_file):
        add_kegg(kegg_file)

    return ret

def fix_img_gff_errors(gff_file):
    ''' GFF files in the IMG directory can contain errors. This fixes them and returns a file handle to the fixed file

        transforms literal semicolon (;) characters in field 9 to their proper percent encoding (%3B)
        fixes the CRISPR lines
    '''
    new_gff = TemporaryFile()
    with open(gff_file) as fp:
        for ln, line in enumerate(fp):
            if line[0] == '#':
                new_gff.write(line)
            else:
                fields = line.split('\t')
                if fields[2] == 'CRISPR':
                    fields[5] = '.'
                    fields[6] = '?'
                    fields[7] = '.'
                else:
                    attributes = fields[8]
                    attribute_kv = attributes.split(';')
                    new_attributes = []
                    for i in range(len(attribute_kv)):
                        if ',' in attribute_kv[i]:
                            attribute_kv[i] = re.sub(',', '%2C', attribute_kv[i])

                        if '=' not in attribute_kv[i]:
                            new_attributes[-1] += '%3B' + attribute_kv[i]
                        else:
                            new_attributes.append(attribute_kv[i])
                    fields[8] = ';'.join(new_attributes)

                new_gff.write('\t'.join(fields))

    new_gff.seek(0)
    return new_gff

def load_img(db, directory, fetch_taxonomy=False, taxid=None):

    dir, bas = path_split(directory)
    if not bas:
        dir, bas = path_split(dir)

    fasta_file = path_join(directory, bas + '.fna')
    gff_file = path_join(directory, bas + '.gff')
    kegg_file = path_join(directory, bas + '.ko.tab.txt')
    xref_file = path_join(directory, bas + '.xref.tab.txt')

    gff_file = fix_img_gff_errors(gff_file)

    xref_dict = None
    if isfile(kegg_file):
        xref_dict = add_gene_dbxref(xref_file, kegg_file)

    with open(fasta_file) as seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

    saved = []
    for rec in GFF.parse(gff_file, seq_dict ):
        for i in range(len(rec.features)):
            feat = rec.features[i]
            try:
                dbxrefs = xref_dict[feat.qualifiers['ID'][0]]
                try:
                    rec.features[i].qualifiers['db_xref'].extend(dbxrefs)
                except:
                    rec.features[i].qualifiers['db_xref'] = dbxrefs
            except KeyError:
                # this gene has no xrefs
                pass

            #print rec.features[i]

        saved.append(add_sequence_dbxref(rec, taxid, bas))

    db.load(saved, fetch_NCBI_taxonomy=fetch_taxonomy)


def main(args):

    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)
    if args.database_name not in server.keys():
        server.new_database(args.database_name)

    db = server[args.database_name]
    try:
        load_img(db, args.directory, args.tax_lookup, args.taxid)
        server.adaptor.commit()
    except:
        server.adaptor.rollback()
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host')
    parser.add_argument('-u', '--user', help='database user name')
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to', default='localhost')
    parser.add_argument('directory', help='directory containing genome information downloaded from IMG. This makes some assumptions on how files and folders are named:'
            ' The directory is assumed to be named with the IMG genome id (e.g. 642555106) and the files within it should have that same ID as their prefix with .gff and .fna'
            ' as file extensions. This will also load in the KEGG annotations available for the genome into the database')
    parser.add_argument('-t', '--lookup-taxonomy', dest='tax_lookup', help='access taxonomy information on NCBI servers', action="store_true", default=False)
    parser.add_argument('-T', '--taxid', help='supply a ncbi taxonomy id that will be applied to all sequences in the file', default=None)
    args = parser.parse_args()
    main(args)

