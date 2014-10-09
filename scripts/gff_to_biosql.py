#!/usr/bin/env python
"""Load a fasta file of sequences and associated GFF file into BioSQL.

You will need to adjust the database parameters and have a BioSQL database set
up. See:

http://biopython.org/wiki/BioSQL

Depending on the size of the sequences being loaded, you may also get errors on
loading very large chromosome sequences. Updating these options can help:

    set global max_allowed_packet=1000000000;
    set global net_buffer_length=1000000;

Usage:
    gff_to_biosql.py <fasta file> <gff file>
"""
from __future__ import with_statement
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO

from BCBio.GFF import GFFParser

def main(seq_file, gff_file):
    # -- To be customized
    # You need to update these parameters to point to your local database
    # XXX demo example could be swapped to use SQLite when that is integrated
    db_name = "orphan.db"
    biodb_name = 'metagenomic_database'

    print "Parsing FASTA sequence file..."
    with open(seq_file) as seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

    print "Parsing GFF data file..."
    parser = GFFParser()
    recs = parser.parse(gff_file, seq_dict )#, limit_info=limit_info)
    for r in recs:
        print r.features[0]
    #print "Writing to BioSQL database..."
    #server = BioSeqDatabase.open_database(driver="sqlite3",db=db_name)
    #try:
    #    if biodb_name not in server.keys():
    #        server.new_database(biodb_name)
    #    else:
    #        server.remove_database(biodb_name)
    #        server.adaptor.commit()
    #        server.new_database(biodb_name)
    #    db = server[biodb_name]
    #    db.load(recs)
    #    server.adaptor.commit()
    #except:
    #    server.adaptor.rollback()
    #    raise

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print __doc__
        sys.exit()
    main(sys.argv[1], sys.argv[2])
