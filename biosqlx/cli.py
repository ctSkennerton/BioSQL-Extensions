# -*- coding: utf-8 -*-

"""Console script for biosqlx."""
import sys
import csv
import click
import os
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
from biosqlx.biosqlx import CustomDBLoader
from dotenv import load_dotenv, find_dotenv

# global server object to be initialized by the main command
# and utilized by all of the subcommands
server = None


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


@click.group()
@click.version_option()
@click.option('-d', '--database', help='name of premade biosql database', default='biosqldb', show_default=True)
@click.option('-r', '--driver', help='Python database driver to use (must be installed separately)',
        type=click.Choice(["MySQLdb", "psycopg2", "sqlite3"]), default='psycopg2', show_default=True)
@click.option('-p', '--port', help='post to connect to on the host',
        default=5432, show_default=True)
@click.option('-u', '--user', help='database user name',
        default=None, show_default=True)
@click.option('-P','--password', help='database password for user')
@click.option('-H', '--host', help='host to connect to',
        default='localhost', show_default=True)
def main(database, driver, port, user, password, host):
    """Console script for biosqlx."""
    global server

    # look in the users home directory for a config file
    dotenv_path = os.path.join(os.path.expanduser("~"), '.biosqlx.cfg')
    # load up the entries as environment variables
    load_dotenv(dotenv_path)

    if not user:
        user = os.environ.get("BIOSQLX_USER", getuser())

    if not password:
        password = os.environ.get("BIOSQLX_PASSWORD")
        if not password:
            password = getpass("Please enter the password for user " + \
                                user + " on database " + database + ": ")

    server = BioSeqDatabase.open_database(driver=driver,
                db=database, user=user, host=host, passwd=password)
    return 0


@main.group()
def add():
    '''Add information to the database'''

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
@click.option('-D', '--database-name', help='limit the extracted sequences from this namespace', default=None)
@click.argument('new_taxons', nargs=-1 )
def sequence(fasta, gff, genbank, lookup_taxonomy, taxid, database_name, new_taxons):
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
    def _add_taxid(inIter, taxid):
        inIter.annotations['ncbi_taxid'] = taxid
        return inIter

    def _load_gff(db, gff_file, fasta_file, fetch_taxonomy=False, taxid=None):
        from BCBio import GFF
        with open(fasta_file) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

        saved = []
        for rec in GFF.parse(gff_file, seq_dict ):
            saved.append(_add_taxid(rec, taxid))

        db.load(saved, fetch_NCBI_taxonomy=fetch_taxonomy)

    def _load_genbank(db, genbank_file, fetch_taxonomy=False, taxid=None):
        with open(genbank_file) as fp:
            saved = []
            for rec in SeqIO.parse(fp, 'genbank' ):
                rec = _add_taxid(rec, taxid)
                saved.append(rec)
            db.load(saved, fetch_NCBI_taxonomy=fetch_taxonomy)

    if database_name not in server.keys():
        server.new_database(database_name)

    db = server[database_name]

    try:
        if gff is not None and fasta is not None:
            _load_gff(db, gff, fasta, lookup_taxonomy, taxid)
        elif genbank is not None:
            _load_genbank(db, genbank, lookup_taxonomy, taxid)
    except e:
        click.echo(e)
        server.adaptor.rollback()
        click.echo("problem loading new records into database",
                file=sys.stderr)
        sys.exit(1)

    if new_taxons:
        taxon_tree = TaxonTree(server.adaptor)
        nodes = taxon_tree.find_elements(ncbi_taxon_id=taxid)
        if len(nodes) == 0:
            # parent doesn't exist, error
            click.echo("The supplied NCBI taxonomy id via -T "
                    "does not appear to be valid. Cannot find "
                    "it in the database. Use -t in conjunction "
                    "to query NCBI servers and add this into "
                    "the database.", file=sys.stderr)
            sys.exit(1)
        elif len(nodes) > 1:
            # name insn't unique, error
            click.echo("There is more than one taxon with the given "
                    "identifier", file=sys.stderr)
            for node in nodes:
                ckick.echo(node, file=sys.stderr)
            sys.exit(1)
        else:
            parent_node = nodes[0]

        new_taxons = map(lambda x: x.split(":"), new_taxons)
        for taxname, taxrank in new_taxons:
            nodes = taxon_tree.find_elements(name=taxname)
            if len(nodes) == 0:
                # this guy doesn't exist yet, add him in
                # and make it the new parent for the next round
                parent_node = taxon_tree.add(taxname, 'scientific name',
                        rank=taxrank, parent=parent_node)
            elif len(nodes) > 1:
                # name insn't unique, error
                pass
            else:
                # this guy already exists
                # so we shouldn't add him in again
                parent_node = nodes[0]


        if fasta is not None:
            gen = SeqIO.parse(fasta, 'fasta')
        elif genbank is not None:
            gen = SeqIO.parse(genbank, 'genbank')

        for rec in gen:
            server.adaptor.execute('UPDATE bioentry SET taxon_id = %s WHERE bioentry_id = %s',
                    (parent_node._id, db.adaptor.fetch_seqid_by_display_id(db.dbid, rec.name)))

    server.commit()


@main.group()
def modify():
    '''Modify existing data in the database'''

@modify.command()
@click.option('-i', '--input', help='provide text file, tab delimited, '
        'where the first column is the name of the sequence feature and '
        'the following columns are the values of the annotation that you '
        'want to add. The first line must be a header line, which will '
        'be used as the name of the qualifier for the seqfeature.')
@click.option('-g', '--gff', help='provide a gff3 formatted file whose '
        'attributes will be added to existing sequence features. Only the '
        'information in the last column of the gff file will be utilized '
        'so you must make sure that either the ID or locus_tag qualifiers '
        'are present in the gff file. If both are present then ID will be '
        'preferred over locus_tag. If neither are present then the record '
        'will be skipped. Make sure that the ID or locus_tag are unique '
        '(and present) in the database otherwise the attributes will not '
        'be loaded.')
@click.option('-s', '--seqfeature', help='The key column of the input '
        'file is the seqfeature id used by the database. Does not apply '
        'when using a gff file as input', is_flag=True, default=False)
@click.option('--replace', help='replace any existing annotations for the '
        'given qualifiers', is_flag=True, default=False)
@click.option('--key', help='name of the column that contains a unique '
        'identifier for the seqfeature. e.g. locus_tag', required=True)
def annotation(infile, gff, seqfeature, replace, key):
    '''Change of add in annotations to existing sequences'''

    def _parse_input(infile, key):
        mapping = {}
        with open(infile) as fp:
            reader = csv.DictReader(fp, delimiter="\t")
            for row in reader:
                mapping[(key, row[key])] = {}
                for k, v in row.items():
                    if k != key:
                        try:
                            mapping[(key, row[key])][k].append(v)
                        except KeyError:
                            mapping[(key, row[key])][k] = [v]
        return mapping

    def _parse_gff(infile):
        from BCBio import GFF
        mapping = {}
        with open(infile) as fp:
            for rec in GFF.parse(fp):
                for feature in rec.features:
                    try:
                        protein_id = feature.qualifiers['ID'][0]
                        mapping[('ID', protein_id)] = feature.qualifiers
                    except KeyError:
                        try:
                            protein_id = feature.qualifiers['locus_tag'][0]
                            mapping[('locus_tag', protein_id)] = feature.qualifiers
                        except KeyError:
                            pass
        return mapping


    db = server[list(server.keys())[0]]

    if infile is not None:
        mapping = _parse_input(infile, key)
    else:
        mapping = _parse_gff(gff)

    # this may be a little tricky depending on how the database is set up
    # since a bioentry is equivelent to a genbank file but genbank files could
    # be created from a whole chromosome or from an individual protein.
    # If it is from a single protein then the protein ID will be the bioentry_id
    # but if it is from a whole genome then it will be a seqfeature_id
    db_loader = CustomDBLoader(db.adaptor, db.dbid, False)

    for (term_name, protein), values in mapping.items():
        # Start by looking for bioentries that have the name
        if not seqfeature:
            seqfeature_id = get_seqfeature_id_from_qv(db, term_name, protein)
        else:
            seqfeature_id = int(protein)
        # now add in our qualifier and value onto that seqfeature
        if replace:
            for qualifier in values.keys():
                if qualifier == 'db_xref':
                    click.echo('Cannot remove any current db_xref, '
                               'you must do this manually for seqfeature {}'.format(seqfeature_id),
                                                                                    file=sys.stderr)
                else:
                    db.adaptor.execute("DELETE FROM seqfeature_qualifier_value WHERE term_id = \
                                (SELECT term_id FROM term WHERE ontology_id = \
                                    (SELECT ontology_id FROM ontology WHERE name = 'Annotation Tags')\
                                AND name = %s)\
                            AND seqfeature_id = %s", (qualifier, seqfeature_id))

        try:
            db_loader._load_seqfeature_qualifiers(values, seqfeature_id)
        except:
            click.echo("Fatal Error: failed to load {} with values {}".format(protein, values), file=sys.stderr)
            click.echo("Check the input file for possible errors", file=sys.stderr)
            sys.exit(1)
    server.commit()

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
@click.option('-T', '--feature-type', help='restrict the results to feature '
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
@click.option('-D', '--database-name', help='limit the extracted sequences from this namespace', default=None)
def sequence(output_format, split_species, feature_type, fuzzy, qualifier, value, taxonomy, database_name):
    '''Extract information about sequences from the database'''
    if qualifier is None and value is None and taxonomy is None:
        click.echo("please provide at least -t (extract by taxonomy) or both -q & -v (qualifier and value)")
        sys.exit(1)

    if feature_type:
        feature_type = list(feature_type)
    elif output_format == 'feat-prot':
        feature_type = ['CDS']
    elif output_format == 'feat-nucl':
        feature_type = ['CDS', 'rRNA', 'tRNA']

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
            for dbname, dbid, seqfeatureid, taxid in seqfeature_bioentries:
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
            if database_name is not None and dbname != database_name:
                continue
            dbids[(dbid, dbname)] = taxon_id

    output_files = {}

    # selecting a qualifier and value is going to be a more specific search, so start there
    if qualifier and value:
        seqfeature_ids = get_seqfeature_ids_from_qv(server, qualifier, value,
                                                    fuzzy=fuzzy, namespace=database_name)

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
