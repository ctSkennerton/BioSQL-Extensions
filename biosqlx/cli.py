# -*- coding: utf-8 -*-

"""Console script for biosqlx."""
import sys
import click
from getpass import getuser, getpass
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from biosqlx.util import print_feature_qv_csv, extract_feature_sql, get_bioentries_from_taxonomy, get_seqfeature_ids_for_bioseqs
from biosqlx.taxon_tree import TaxonTree


# global server object to be initialized by the main command
# and utilized by all of the subcommands
server = None

@click.group()
@click.version_option()
@click.option('-d', '--database', help='name of premade biosql database', default='biosqldb', show_default=True)
@click.option('-D', '--database-name', help='limit the extracted sequences from this namespace', default=None)
@click.option('-r', '--driver', help='Python database driver to use (must be installed separately)',
        type=click.Choice(["MySQLdb", "psycopg2", "sqlite3"]), default='psycopg2', show_default=True)
@click.option('-p', '--port', help='post to connect to on the host',
        default=5432, show_default=True)
@click.option('-u', '--user', help='database user name',
        default=getuser(), show_default=True)
@click.option('-P','--password', help='database password for user')
@click.option('-H', '--host', help='host to connect to',
        default='localhost', show_default=True)
def main(database, database_name, driver, port, user, password, host):
    """Console script for biosqlx."""
    global server
    if password is None:
        password = getpass("Please enter the password for user " + \
                user + " on database " + database + ": ")

    server = BioSeqDatabase.open_database(driver=driver,
                db=database, user=user, host=host, passwd=password)
    return 0


@main.group()
def add():
    '''Add information to the database'''
    click.echo("help")

@add.command()
@click.option('-f', '--fasta', help='fasta file to add into the database')
@click.option('-g', '--gff', help='gff file of features to add into the database. '
        'Must be paired with a fasta file')
@click.option('-G', '--genbank', help='genbank file to add into the database')
@click.option('-t', '--lookup-taxonomy', help='Certain files such as genbank files '
        'contain taxonomic information about the seuquence. This option specifies that '
        'this taxonomic information will be queried and downloaded from the NCBI servers',
        is_flag=True, default=False)
@click.option('-T', '--taxid', help='supply a ncbi taxonomy id that will be applied to '
        'all sequences in the file, or if new_taxons are supplied on the command line '
        'this taxonomy ID will be used as the parent taxonomy for the novel lineages. '
        'An error will occur if this taxid is not present in the database and --lookup-taxonomy is false.',
        default=None)
@click.argument('new_taxons', nargs=-1 )
def sequence(fasta, gff, genbank, lookup_taxonomy, taxid, new_taxons):
    '''Add new sequences to the database.

        The New_taxons arguments allow you to specify novel taxonomies
        not currenly in the NCBI database; particularly useful when
        adding in genome bins or novel isolates yet to make it to public
        databases.

        Each taxon specified on the command line should take the
        form of <taxon_name>:<taxon_rank>. Check the taxon table
        in the database for the appropriate values for the taxon_rank,
        however in general they will be names you are familiar with:
        phylum, class, order, family, genus, species
        e.g. ANME-2ab:family ANME-2b:genus ANME-hr1:species

    '''


@main.group()
def modify():
    '''Modify existing data in the database'''

@modify.command()
def annotation():
    '''Change of add in annotations to existing sequences'''

@main.group()
def export():
    '''Extract information from the database'''


@export.command()
@click.option('-o', '--output-format', help='output format of the selected sequences. '
        'Choices: fasta - fasta file of the contigs; gb - genbank file of the sequences; '
        'feat-prot - fasta file containing the translated coding sequences; '
        'feat-nucl - fasta file containing the untranslated coding sequences, '
        'tRNAs and rRNAs; csv - csv file of annotations for the features',
        type=click.Choice(['fasta', 'gb', 'feat-prot', 'feat-nucl', 'csv']), default='fasta')
@click.option('-s', '--split-species', help='when there are multiple species to '
        'be returned, split them into separate files, based on their name, '
        'instead of printing to stdout', is_flag=True, default=False)
@click.option('-t', '--feature-type', help='restrict the results to feature '
        'type e.g. rRNA, tRNA, CDS. This option can be specified multiple '
        'times for multiple types. This option is only used when the ouput format is '
        'feat-nucl or feat-prot, inwhich case the default values will be: CDS, tRNA, rRNA '
        'for feat-nucl; or CDS for feat-prot', default=None, multiple=True)
@click.option('-f', '--fuzzy', help='the value can be a partial match', is_flag=True, default=False)
@click.option('-q', '--qualifier', help='name of the qualifier')
@click.option('-v', '--value', help='value to match on' )
@click.option('-t', '--taxonomy', help='supply a taxonomy name that will be extracted. '
        'If an integer is supplied it will be interpreted as an NCBI '
        'taxonomy id; otherwise it will be interpreted as part of a taxonomy name (e.g. Proteobacteria)')
def sequence(output_format, split_species, feature_type, fuzzy, qualifier, value, taxonomy):
    '''Extract information about sequences from the database'''
    click.echo('me')


    rows = get_bioentries_from_taxonomy(server, taxonomy)
    if len(rows) == 0:
        click.echo("\nThere does not appear to be any sequences associated with\n"
                "the taxonomy provided. If you used a taxonomy name, make sure\n"
                "it is spelled correctly. And remember that it must be the complete name\n"
                "for a particular rank, for example 'Deltaproteo' will match nothing\n"
                "it has to be 'Deltaproteobacteria'.\n"
                "Don't forget to add 'Candidatus ' to the begining of some names\n"
                "or the strain designation for a species. If you used an NCBI taxonomy ID, make\n"
                "sure that it is correct by double checking on the NCBI taxonomy website.", file=sys.stderr)
        sys.exit(1)


    dbids = {}
    for row in rows:
        dbids[(row[0], row[2])] = row[1]
    files = {}
    taxid_to_dbids = {}
    if split_species:
        taxon_file_mapping = {}
        for k, v in dbids.items():
            tname = server.adaptor.execute_and_fetch_col0(
                    "SELECT name from taxon_name where taxon_id = %s and name_class = %s",
                    (v,'scientific name'))[0]
            tname = tname.replace(' ', '_')
            if output_format == 'gb':
                tname += '.gb'
            elif output_format == 'feat-prot':
                tname += '.faa'
            elif output_format == 'csv':
                tname += '.csv'
            else:
                tname += '.fna'
            files[v] = tname
            taxid_to_dbids.setdefault(v, []).append(k)


    if split_species:
        # got to save all of the records before printing them out
        outdata = {}
        for taxid, dbid_list in taxid_to_dbids.items():
            if output_format == 'csv':
                with open(files[taxid], 'w') as fp:
                    print_feature_qv_csv(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbid_list]), fp)
            else:
                for dbid, dbname in dbid_list:
                    db = server[dbname]
                    seq_rec = db[dbid]
                    outdata.setdefault(taxid, []).append(seq_rec)

        for taxid, dbrecs in outdata.items():
            with open(files[taxid], 'w') as fp:
                if 'feat' in output_format:
                    for dbrec in dbrecs:
                        extract_feature(dbrec, output_format, fp)
                elif 'csv' != output_format:
                    SeqIO.write(dbrecs, fp, output_format)

    else:

        if feature_type is not None:
            types = feature_type
        elif output_format == 'feat-prot':
            types = ['CDS']
        elif output_format == 'feat-nucl':
            types = ['CDS', 'rRNA', 'tRNA']

        if output_format == 'feat-prot':
            extract_feature_sql(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]),type=types, translate=True )
        elif output_format == 'feat-nucl':
            extract_feature_sql(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]), type=types)
        elif output_format == 'csv':
            print_feature_qv_csv(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]))
        else:
            for (dbid, dbname), taxid in dbids.items():
                db = server[dbname]
                try:
                    dbrec = db[dbid]
                    SeqIO.write(dbrec, sys.stdout, output_format)
                except KeyError:
                    pass

@export.command()
@click.option('-r', '--root', help='Specify the root of the output tree. '
        'The default is to print all of the organisms in the tree, but by '
        'using this option you can specify a subtree to print', default=None)
def taxonomy(root):
    '''Get information about the organisms present in the database'''
    tree = TaxonTree(server.adaptor)
    if root is not None:
        elements = tree.find_elements(name=root)
        if len(elements) != 1:
            click.echo("The name {} is found more than once or not at all; cannot use as the root node".format(root), file=sys.stderr)
            sys.exit(1)
        else:
            root = elements[0]
    click.echo(tree.pretty_print(root_node=root))


@main.group()
def info():
    '''Get information about the database'''
    click.echo('bleg')


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
