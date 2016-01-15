#!/usr/bin/env python

import argparse
import yaml
from BioSQL import BioSeqDatabase



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=argparse.FileType(),
            help='Input file comes from kegg brite database that has been processed with keg2yml.py')
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-D', '--database-name', help='namespace of the database that you want to add into', dest='database_name', default='metagenomic_database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host')
    parser.add_argument('-u', '--user', help='database user name')
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-o', '--ontology', help='give a name for this ontology', default='KEGG')

    args = parser.parse_args()
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)


    # create an ontology for KEGG if it doesn't already exist
    rows = server.adaptor.execute_and_fetchall('select ontology_id from ontology where name = %s', (args.ontology))
    if len(rows) == 0:
        server.adaptor.execute('insert into ontology (name, definition) values(%s, %s)', (args.ontology, 'The Kegg ortholog ontology'))
        kegg_ontology = server.adaptor.execute_and_fetch_col0('select ontology_id from ontology where name = %s', (args.ontology))[0]
    else:
        kegg_ontology = rows[0][0]

    # Get the predicate term id
    rows = server.adaptor.execute_and_fetchall('select term_id from term where name = "partOf"')
    if len(rows) == 0:
        server.adaptor.execute('insert into term (name, definition) values(%s, %s)', ('partOf', 'Used as a predicate when defining the relationship between terms'))
        predicate_id = server.adaptor.execute_and_fetch_col0('select term_id from term where name = "partOf"')[0]
    else:
        predicate_id = rows[0][0]

    yaml_data = yaml.load(args.infile)
    for l1,d1 in yaml_data.items():
        print l1
        server.adaptor.execute('insert into term (name, ontology_id) values(%s, %s) ', (l1, kegg_ontology))
        l1_id = server.adaptor.execute_and_fetch_col0('select term_id from term where name = %s and ontology_id = %s', (l1, kegg_ontology))[0]

        for l2,d2 in d1.items():
            print '    -',l2
            server.adaptor.execute('insert into term (name, ontology_id) values(%s,%s)', (l2, kegg_ontology))
            l2_id = server.adaptor.execute_and_fetch_col0('select term_id from term where name = %s and ontology_id = %s', (l2, kegg_ontology))[0]
            server.adaptor.execute('insert into term_relationship (subject_term_id, predicate_term_id, object_term_id, ontology_id) values (%s, %s, %s, %s)', (l1_id, predicate_id, l2_id, kegg_ontology))
            server.adaptor.execute('insert into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l1_id, predicate_id, l2_id, kegg_ontology, 1))

            for l3,d3 in d2.items():
                if d3['orthologs'] is None:
                    continue

                print '        -',l3,d3['definition']
                server.adaptor.execute('insert into term (identifier, name, ontology_id) values (%s, %s, %s)', (l3, d3['definition'], kegg_ontology))
                l3_id = server.adaptor.execute_and_fetch_col0('select term_id from term where name = %s and ontology_id = %s', (l3, kegg_ontology))[0]
                server.adaptor.execute('insert into term_relationship (subject_term_id, predicate_term_id, object_term_id, ontology_id) values (%s, %s, %s, %s)', (l2_id, predicate_id, l3_id, kegg_ontology))
                server.adaptor.execute('insert into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l1_id, predicate_id, l3_id, kegg_ontology, 2))
                server.adaptor.execute('insert into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l2_id, predicate_id, l3_id, kegg_ontology, 1))

                for ortholog in d3['orthologs']:
                    # check if we already have this kegg ortholog
                    rows = server.adaptor.execute_and_fetchall('select term_id from term where identifier = %s', (ortholog,))
                    if len(rows) != 0:
                        ortho_id = rows[0][0]
                    else:
                        server.adaptor.execute('insert into term (identifier, name, ontology_id) values(%s,%s,%s)', (ortholog, ortholog, kegg_ontology))
                        ortho_id = server.adaptor.execute_and_fetch_col0('select term_id from term where identifier = %s and name = %s and ontology_id = %s', (ortholog, ortholog, kegg_ontology) )[0]

                    server.adaptor.execute('insert or ignore into term_relationship (subject_term_id, predicate_term_id, object_term_id, ontology_id) values (%s, %s, %s, %s)', (l3_id, predicate_id, ortho_id, kegg_ontology))
                    server.adaptor.execute('insert or ignore into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l1_id, predicate_id, ortho_id, kegg_ontology, 3))
                    server.adaptor.execute('insert or ignore into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l2_id, predicate_id, ortho_id, kegg_ontology, 2))
                    server.adaptor.execute('insert or ignore into term_path (subject_term_id, predicate_term_id, object_term_id, ontology_id, distance) values (%s,%s,%s,%s,%s)', (l3_id, predicate_id, ortho_id, kegg_ontology, 1))

    server.commit()



