#!/usr/bin/env python
import argparse
import yaml
import sys
from BioSQL import BioSeqDatabase

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host')
    parser.add_argument('-u', '--user', help='database user name')
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-i', '--infile', type=argparse.FileType(), help='yaml file containing information about kegg orthologs')
    args = parser.parse_args()

    data = yaml.load(args.infile)
    #print yaml.dump(yaml.load(args.infile), default_flow_style=False)
    #sys.exit(0)
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)
    dbxref_check = 'select * from dbxref where accession = %s'
    dbxref_insert = 'insert into dbxref (dbname, accession, version) values ("ko", %s, 1)'
    term_dbxref_insert = 'insert into term_dbxref (term_id, dbxref_id) values (%s, %s)'
    term_id_ko = 'select term_id from term where identifier = %s'
    ko_ontology_id_select = "select ontology_id from ontology where name = %s"
    ko_ontology_insert = "insert into ontology (name) values (%s)"
    term_ko_insert = 'insert into term (name, identifier, ontology_id) values (%s, %s, %s)'
    dbxref_qv_insert = 'insert into dbxref_qualifier_value values (%s, %s, %s, %s)'


    kegg_ontology_id = server.adaptor.execute_and_fetch_col0(ko_ontology_id_select, 'KEGG')[0]
    if not kegg_ontology_id:
        server.adaptor.execute(ko_ontology_insert, 'KEGG')
        kegg_ontology_id = server.adaptor.execute_and_fetch_col0(ko_ontology_id_select, 'KEGG')[0]


    for ko, qv in data.items():
        # get the term id for this KO
        rows = server.adaptor.execute_and_fetchall(term_id_ko, (ko,))
        if len(rows) == 0:
            # there are no terms that match this KO, need to add it in
            server.adaptor.execute(term_ko_insert, (ko, ko, kegg_ontology_id))
            ko_term_id = server.adaptor.execute_and_fetchall(term_id_ko, (ko,))
        elif len(rows) == 1:
            ko_term_id = rows[0][0]
        else:
            raise RuntimeError("multiple terms with "+ko)

        rows = server.adaptor.execute_and_fetchall(dbxref_check, (ko,))
        if len(rows) == 1:
            dbxref_id = rows[0][0]
        elif len(rows) > 1:
            raise RuntimeError("multiple dbxrefs for "+ko)
        else:
            # no dbxref
            # now insert the crossref and link it to the term
            server.adaptor.execute(dbxref_insert, (ko,))
            dbxref_id = server.adaptor.execute_and_fetch_col0(dbxref_check, (ko,))[0]

        rows = server.adaptor.execute_and_fetchall('select * from term_dbxref where dbxref_id = %s and term_id = %s', (dbxref_id, ko_term_id))
        if len(rows) == 0:
            try:
                server.adaptor.execute(term_dbxref_insert, (ko_term_id, dbxref_id))
            except sqlite3.IntegrityError, e:
                print ko_term_id, dbxref_id
                raise e


        # add in the qualifiers for the cross reference
        for term_name, value in qv.items():
            term_id = server.adaptor.execute_and_fetch_col0('select term_id from term where name = %s', (term_name,))[0]
            rank = 0

            if isinstance(value, list):
                for i in value:
                    rows = server.adaptor.execute_and_fetchall('select * from dbxref_qualifier_value where dbxref_id = %s and term_id = %s and value = %s', (dbxref_id, term_id, i))
                    if len(rows) == 0:
                        server.adaptor.execute(dbxref_qv_insert, (dbxref_id, term_id, rank, i))
                    rank += 1
            else:
                rows = server.adaptor.execute_and_fetchall('select * from dbxref_qualifier_value where dbxref_id = %s and term_id = %s and value = %s', (dbxref_id, term_id, i))
                if len(rows) == 0:
                    server.adaptor.execute(dbxref_qv_insert, (dbxref_id, term_id, rank, value))

    server.commit()


