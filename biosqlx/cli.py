# -*- coding: utf-8 -*-

"""Console script for biosqlx."""
import sys
import csv
import click
from getpass import getuser, getpass
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from biosqlx.util import print_feature_qv_csv, \
        extract_feature_sql, \
        get_bioentries_from_taxonomy, \
        get_seqfeature_ids_for_bioseqs, \
        get_seqfeature_ids_from_qv, \
        get_bioseqid_for_seqfeature
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
    if qualifier is None and value is None and taxonomy is None:
        click.echo("please provide at least -t (extract by taxonomy) or both -q & -v (qualifier and value)")
        sys,exit(1)

    if feature_type:
        feature_type = list(feature_type)
    elif output_format == 'feat-prot':
        feature_type = ['CDS']
    elif output_format == 'feat-nucl':
        feature_type = ['CDS', 'rRNA', 'tRNA']

    def _check_tax(server, taxonomy):
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
        else:
            return rows


    def _make_file_mapping(server, dbids):
        files = {}
        taxid_to_dbids = {}
        for k, v in dbids.items():
            if v is None:
                # special case for no taxonomy
                tname = 'Unassigned'
            else:
                try:
                    tname = server.adaptor.execute_and_fetch_col0(
                            "SELECT name from taxon_name where taxon_id = %s and name_class = %s",
                            (v,'scientific name'))[0]
                except IndexError:
                    raise RuntimeError("cannot get the scientific name for {}".format(v))
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
        return files, taxid_to_dbids


    def _choose_output_format(server, sfids, feature_type, output_format, ofile=sys.stdout, bioentries=None):
        if output_format == 'feat-prot':
            extract_feature_sql(server, sfids, type=feature_type, translate=True, file=ofile )
        elif output_format == 'feat-nucl':
            extract_feature_sql(server, sfids, type=feature_type, file=ofile)
        elif output_format == 'csv':
            print_feature_qv_csv(server, sfids, outfile=ofile)
        else:
            # the user wants the fasta or genbank options, in which case the
            # whole biosequence will be returned
            seqfeature_bioentries = get_bioseqid_for_seqfeature(server, sfids)
            printed_bioseqs = set()
            for dbname, dbid, seqfeatureid in seqfeature_bioentries:
                if dbid not in printed_bioseqs:
                    db = server[dbname]
                    try:
                        dbrec = db[dbid]
                        SeqIO.write(dbrec, ofile, output_format)
                    except KeyError:
                        pass
                    finally:
                        printed_bioseqs.add(dbid)

    dbids = {}
    if taxonomy:
        rows = _check_tax(server, taxonomy)
        for dbid, taxon_id, dbname in rows:
            dbids[(dbid, dbname)] = taxon_id

    output_files = {}

    # selecting a qualifier and value is going to be a more specific search, so start there
    if qualifier and value:
        seqfeature_ids = get_seqfeature_ids_from_qv(server, qualifier, value, fuzzy=fuzzy)

        # if the output_format is just features and is all going to one output file then we
        # can move on here to printing
        if taxonomy is None:
            if split_species:
                seqfeature_bioentries = get_bioseqid_for_seqfeature(server, seqfeature_ids)
                dbid_to_seqfeature_id = {}
                filtered_tax = {}
                for dbname, dbid, seqfeatureid, taxon_id in seqfeature_bioentries:
                    filtered_tax[(dbid, dbname)] = taxon_id
                    try:
                        dbid_to_seqfeature_id[dbid].append(seqfeatureid)
                    except KeyError:
                        dbid_to_seqfeature_id[dbid] = [seqfeatureid]

                files, taxid_to_dbid = _make_file_mapping(server, filtered_tax)
                taxid_to_seqfeature = {}
                for taxid, dbid_list in taxid_to_dbid.items():
                    taxid_seqfeatures = []
                    for dbid, dbname in dbid_list:
                        taxid_seqfeatures.extend(dbid_to_seqfeature_id[dbid])

                    with open(files[taxid], 'w') as fp:
                        _choose_output_format(server, taxid_seqfeatures,
                                              feature_type, output_format,
                                              ofile=fp, bioentries=None)
            else:
                _choose_output_format(server, seqfeature_ids, feature_type, output_format)
                #seqfeature_bioentries = get_bioseqid_for_seqfeature(server, seqfeature_ids)
                #for dbname, dbid, seqfeatureid in seqfeature_bioentries:
                #    taxid = dbids[(dbid, dbname)]
                #    outfile_name = _make_file_mapping(server, taxid)
                #    if outfile_name not in output_files:
                #        output_files[outfile_name] = open(outfile_name, 'w')

                #    _choose_output_format(server, seqfeature_ids, feature_type, output_format, output_files[outfile_name])
        else:
            # filter the seqfeatures based on the wanted taxonomy
            seqfeature_bioentries = get_bioseqid_for_seqfeature(server, seqfeature_ids)
            final_seqfeatures = []
            dbid_to_seqfeature_id = {}
            filtered_tax = {}
            for dbname, dbid, seqfeatureid, taxon_id in seqfeature_bioentries:
                if (dbid, dbname) in dbids:
                    filtered_tax[(dbid, dbname)] = dbids[(dbid, dbname)]
                    final_seqfeatures.append(seqfeatureid)
                    try:
                        dbid_to_seqfeature_id[dbid].append(seqfeatureid)
                    except KeyError:
                        dbid_to_seqfeature_id[dbid] = [seqfeatureid]

            if split_species:
                # get a mapping of the taxonomy
                files, taxid_to_dbid = _make_file_mapping(server, filtered_tax)
                taxid_to_seqfeature = {}
                for taxid, dbid_list in taxid_to_dbid.items():
                    taxid_seqfeatures = []
                    for dbid, dbname in dbid_list:
                        taxid_seqfeatures.extend(dbid_to_seqfeature_id[dbid])

                    with open(files[taxid], 'w') as fp:
                        _choose_output_format(server, taxid_seqfeatures,
                                              feature_type, output_format,
                                              ofile=fp, bioentries=None)
            else:
                _choose_output_format(server, final_seqfeatures, feature_type, output_format)
    else:
        # no qualifier and value
        if split_species:
            dbid_to_seqfeature_id = {}
            for dbid, dbname in dbids.keys():
                dbid_to_seqfeature_id[dbid] = get_seqfeature_ids_for_bioseqs(server, [dbid])

            # get a mapping of the taxonomy
            files, taxid_to_dbid = _make_file_mapping(server, dbids)
            taxid_to_seqfeature = {}
            for taxid, dbid_list in taxid_to_dbid.items():
                taxid_seqfeatures = []
                for dbid, dbname in dbid_list:
                    taxid_seqfeatures.extend(dbid_to_seqfeature_id[dbid])

                with open(files[taxid], 'w') as fp:
                    _choose_output_format(server, taxid_seqfeatures,
                                          feature_type, output_format,
                                          ofile=fp, bioentries=None)
        else:
            final_seqfeatures = get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()])
            _choose_output_format(server, final_seqfeatures, feature_type, output_format)


@export.command()
@click.option('-r', '--root', help='Specify the root of the output tree. '
        'The default is to print all of the organisms in the tree, but by '
        'using this option you can specify a subtree to print', default=None)
@click.option('-o', '--output-format', type=click.Choice(['tree', 'lineage']),
        help='choose differnt output formats for the data. The tree format '
        'will print a nice hierarchical representation; lineage will print a '
        'semicolon (;) separated list of the taxonomy', default='tree')
def taxonomy(root, output_format):
    '''Get information about the organisms present in the database'''
    tree = TaxonTree(server.adaptor)
    if root is not None:
        elements = tree.find_elements(name=root)
        if len(elements) != 1:
            click.echo("The name {} is found more than once or not at all; cannot use as the root node".format(root), file=sys.stderr)
            sys.exit(1)
        else:
            root = elements[0]
    if output_format == 'tree':
        click.echo(tree.pretty_print(root_node=root))
    elif output_format == 'lineage':
        writer = csv.writer(sys.stdout)
        for row in tree.lineage(root_node=root):
            writer.writerow(row)



@main.group()
def info():
    '''Get information about the database'''
    click.echo('bleg')


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
