#!/usr/bin/env python
import sys
import argparse
from BioSQL import BioSeqDatabase
from Bio import SeqIO

def add_taxid(inIter, taxid):
    inIter.annotations['ncbi_taxid'] = taxid
    yield inIter

def load_gff(db, gff_file, fasta_file, fetch_taxonomy=False, taxid=None):
    from BCBio.GFF import GFFParser
    with open(fasta_file) as seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

    parser = GFFParser()
    recs = parser.parse(gff_file, seq_dict )#, limit_info=limit_info)
    db.load(add_taxid(recs, taxid), fetch_NCBI_taxonomy=fetch_taxonomy)

def load_genbank(db, genbank_file, fetch_taxonomy=False, taxid=None):
    with open(genbank_file) as fp:
        db.load(add_taxid(SeqIO.parse(genbank_file, 'genbank'), taxid), fetch_NCBI_taxonomy=fetch_taxonomy)


def main(args):
    server = BioSeqDatabase.open_database(driver="sqlite3",db=args.database)
    if args.database_name not in server.keys():
        server.new_database(args.database_name)

    db = server[args.database_name]
    try:
        if args.gff is not None and args.fasta is not None:
            load_gff(db, args.gff, args.fasta, args.tax_lookup)
            server.adaptor.commit()
        elif args.genbank is not None:
            load_genbank(db, args.genbank, args.tax_lookup)
            server.adaptor.commit()
    except:
        server.adaptor.rollback()
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-f', '--fasta', help='fasta file to add into the database')
    parser.add_argument('-g', '--gff', help='gff file of reatures to add into the database. Must be paired with a fasta file')
    parser.add_argument('-G', '--genbank', help='genbank file to add into the database')
    parser.add_argument('-t', '--lookup-taxonomy', dest='tax_lookup', help='access taxonomy information on NCBI servers', action="store_true", default=False)
    parser.add_argument('-T', '--taxid', help='supply a ncbi taxonomy id that will be applied to all sequences in the file')
    args = parser.parse_args()
    main(args)

