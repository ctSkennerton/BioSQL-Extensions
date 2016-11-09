#!/usr/bin/env python
import sys
import argparse
from getpass import getpass, getuser
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'c.skennerton@gmail.com'

def add_taxid(inIter, taxid):
    inIter.annotations['ncbi_taxid'] = taxid
    return inIter

def load_gff(db, gff_file, fasta_file, fetch_taxonomy=False, taxid=None):
    from BCBio import GFF
    with open(fasta_file) as seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

    saved = []
    for rec in GFF.parse(gff_file, seq_dict ):
        saved.append(add_taxid(rec, taxid))

    db.load(saved, fetch_NCBI_taxonomy=fetch_taxonomy)

def load_genbank(db, genbank_file, fetch_taxonomy=False, taxid=None):
    with open(genbank_file) as fp:
        saved = []
        for rec in SeqIO.parse(fp, 'genbank' ):
            rec = add_taxid(rec, taxid)
            saved.append(rec)
        db.load(saved, fetch_NCBI_taxonomy=fetch_taxonomy)


def update_left_right_taxon_values(server, left_value):
    """ update the left and right values in the table
    """
    if not left_value:
        return
    # Due to the UNIQUE constraint on the left and right values in the taxon
    # table we cannot simply update them through an SQL statement as we risk
    # colliding values. Instead we must select all of the rows that we want to
    # update, modify the values in python and then update the rows
    #self.adaptor.execute("UPDATE taxon SET right_value = right_value + 2 WHERE right_value >= %s", (left_value,))
    #self.adaptor.execute("UPDATE taxon SET left_value = left_value + 2 WHERE left_value > %s", (left_value,))

    rows = server.adaptor.execute_and_fetchall(
            "SELECT left_value, right_value, taxon_id FROM taxon WHERE right_value >= %s or left_value > %s",
            #"SELECT left_value, right_value, taxon_id FROM taxon",
            (left_value, left_value)
            )

    right_rows = []
    left_rows = []
    for i in range(len(rows)):
        new_right = rows[i][1]
        new_left = rows[i][0]
        if new_right >= left_value:
            new_right += 2

        if new_left > left_value:
            new_left += 2
        right_rows.append((new_right, rows[i][2]))
        left_rows.append((new_left, rows[i][2]))


    # sort the rows based on the value from largest to smallest
    # should ensure no overlaps
    right_rows = sorted(right_rows, key=lambda x: x[0], reverse=True)
    left_rows = sorted(left_rows, key=lambda x: x[0], reverse=True)

    try:
        server.adaptor.executemany("UPDATE taxon SET left_value = %s WHERE taxon_id = %s", left_rows)
        server.adaptor.executemany("UPDATE taxon SET right_value = %s WHERE taxon_id = %s", right_rows)
    except:
        #print(orig_rows)
        #print(rows2)
        #print(rows)

        #for r in self.adaptor.execute_and_fetchall("select * from taxon"):
        #    print(r)

        raise

def insert_taxon_rank(server, parent_taxon_id, parent_left_value, parent_right_value, node_name, node_rank):
    # make sure that the name doesn't already exist
    rows = server.adaptor.execute_and_fetch_col0('select taxon_id from taxon_name where name = %s', (node_name,))
    if rows:
        return server.adaptor.execute_and_fetchall(
                'select taxon_id, left_value, right_value from taxon where taxon_id = %s ',
                (rows[0], ))[0]
    else:
        left_value = parent_right_value
        right_value = parent_right_value + 1
        update_left_right_taxon_values(server, left_value)
        server.adaptor.execute(
                'insert into taxon (parent_taxon_id, node_rank, left_value, right_value)'
                ' values (%s, %s, %s, %s)',
                (parent_taxon_id, node_rank, left_value, right_value))

        taxon_id = server.adaptor.execute_and_fetch_col0(
                'select taxon_id from taxon where parent_taxon_id = %s '
                'and left_value = %s and right_value = %s',
                (parent_taxon_id, left_value, right_value))[0]

        server.adaptor.execute(
            "INSERT INTO taxon_name(taxon_id, name, name_class)"
            " VALUES (%s, %s, 'scientific name')", (taxon_id,
                                                    node_name[:255]))

        return taxon_id, left_value, right_value

def add_new_taxonomy(server, new_taxons, parent_ncbi_tax_id):
    new_taxons = map(lambda x: x.split(":"), new_taxons)
    parent_taxid = None
    parent_left_value = None
    parent_right_value = None
    if parent_ncbi_tax_id:
        parent_taxid, parent_left_value, parent_right_value = \
                server.adaptor.execute_one(
                        'select taxon_id, left_value, right_value '
                        'from taxon where ncbi_taxon_id = %s',
                        (parent_ncbi_tax_id,))

    for tax_name, tax_rank in new_taxons:
        parent_taxid, parent_left_value, parent_right_value = insert_taxon_rank(server,
                parent_taxid,
                parent_left_value,
                parent_right_value,
                tax_name,
                tax_rank)

    return parent_taxid

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)
    if args.database_name not in server.keys():
        server.new_database(args.database_name)

    db = server[args.database_name]

    try:
        if args.gff is not None and args.fasta is not None:
            load_gff(db, args.gff, args.fasta, args.tax_lookup, args.taxid)
        elif args.genbank is not None:
            load_genbank(db, args.genbank, args.tax_lookup, args.taxid)
    except:
        server.adaptor.rollback()
        raise

    if args.new_taxons:
        taxon_id = add_new_taxonomy(server, args.new_taxons, args.taxid)

        if args.fasta is not None:
            gen = SeqIO.parse(args.fasta, 'fasta')
        elif args.genbank is not None:
            gen = SeqIO.parse(args.genbank, 'genbank')

        for rec in gen:
            server.adaptor.execute('update bioentry set taxon_id = %s where bioentry_id = %s',(taxon_id, db.adaptor.fetch_seqid_by_display_id(db.dbid, rec.name)))

    server.commit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host',
            default=5432)
    parser.add_argument('-u', '--user', help='database user name',
            default=getuser())
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to',
            default='localhost')
    parser.add_argument('-f', '--fasta', help='fasta file to add into the database')
    parser.add_argument('-g', '--gff', help='gff file of reatures to add into the database. Must be paired with a fasta file')
    parser.add_argument('-G', '--genbank', help='genbank file to add into the database')
    parser.add_argument('-t', '--lookup-taxonomy', dest='tax_lookup', help='access taxonomy information on NCBI servers', action="store_true", default=False)
    parser.add_argument('-T', '--taxid', help='supply a ncbi taxonomy id that will be applied to all sequences in the file, or if new_taxons are supplied on the command line this taxonomy ID will be used as the parent taxonomy for the novel lineages', default=None)
    parser.add_argument('new_taxons', nargs="*", help='specify novel taxonomies not currenly in the NCBI database. each taxon specified on the command line should take the form of <taxon_name>:<taxon_rank>. Check the taxon table in the database for the appropriate values for the taxon_rank. e.g. ANME-2ab:family ANME-2b:genus ANME-hr1:species')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

