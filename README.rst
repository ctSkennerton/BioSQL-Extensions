.. contents::

Introduction
============

This is a repository of sql tables and views and scripts that could be
added onto an existing BioSQL database. There are a number of scripts
that should hopefully work on any BioSQL database and have been
developed for easy access to various information.

In addition there tables that were made solely for sqlite3 and are
primarily related to adding functionality for performing proteomics
experiments. I have not tried to write the proteomics tables in a
porable way, they were designed to model the data created by the SIPROS
program.

Command line
============

There is a main script ``biosqlx`` that serves as the main extry point
for quickly accessing information in the database. There are a number of
subcommands that allow you to add, modify or export information from the
database. For example ``biosqlx add sequence`` allows you to add new
sequence data into the database, while ``biosqlx export sequence`` lets
you extract data out.

connecting to the database
--------------------------

The top level command ``biosqlx`` contains a few options that let you
change how to connect to the database.

-  ``-d``: This is the name of the database to connect to
-  ``-r``: This is the database driver to use, according to the biosql
   documentation
-  ``-u``: The user name to use to connect to the database.
-  ``-P``: Password for the user to login to the database
-  ``-H``: The host name that the database can be found on
-  ``-p``: The port on the host that the database is found on

Config file
-----------

It's onerous to keep adding in the database user and password for every
command. To alleviate this, ``biosqlx`` can use a configuration file to
store sensitive data like the user name and password. Create a file in
your home directory called ".biosqlx.cfg" and add in the following lines
replacing the template variables:

::

    BIOSQLX_USER=<your_username>
    BIOSQLX_PASSWORD=<your_password>

now change the permisions of that file so only you can read or write to
that file

::

    chmod 600 ~/.biosqlx.cfg


``biosqlx export``
------------------


``biosqlx export sequence``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subcommand allows you to slice and dice the sequence data using
multiple filtering metrics, including the assigned taxonomy or the
sequence feature annotations.

Output formats
^^^^^^^^^^^^^^

The currently supported output formats are: 

1. gb (Genbank), outputs a genbank formatted file for
   the sequence records, these are likely contigs or whole chromosomes.
2. fasta, outputs the sequence information in fasta format. 
3. feat-prot, outputs the translated sequence features (i.e. proteins/CDS)
4. feat-nucl, outputs sequence features untranslated (i.e. proteins/CDS, tRNA, rRNA)
5. csv, output the annotations for all of the features.
6. gff, output annotations as a gff3 formatted file. Just like the genbank
   or fasta output formats, this will print the whole contig.

Filtering on taxonomy
^^^^^^^^^^^^^^^^^^^^^

The input is either an NCBI taxonomy ID or the complete taxonomic name
given with ``-t`` or ``--taxonomy`` on the command line. below are some
examples to show how it works

Extract translated gene sequences from all of the Deltaproteobacteria in
the database

::

    biosqlx export sequence -o feat-prot -t Deltaproteobacteria

Extract all of the sequences from the database that match the NCBI
taxonomy ID, 872. In this case that refers to the
`Desulfovibrio <http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=872&lvl=3&lin=f&keep=1&srchmode=1&unlock>`__.

::

    biosqlx export sequence -o feat-prot -t 872

When using the taxonomic name has the limitation that it must be a
complete name so for example if you wanted sequences from a species then
you must give both the genus and species name. In the below example, if
you just type 'vulgaris', no results will be returned.

::

    biosqlx export sequence -o feat-prot 'Desulfovibrio vulgaris'

Note the single quotes around the organism name! This is important so
that the name is properly input into the script. Any names that are
multiple words need to have the quotes surrounding them.

Filtering using qualifiers and values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You must specify both ``-q`` and ``-v`` on the command line. This is a
strict evaluation of equality so everything must be spelled correctly
etc. Unfortunately this also means that you cannot do other operations
like asking for all features that are greater than or less than some
value. Qualifiers are things like "gene", "product", "db\_xref" and the
associated values would be like "omcX", "cytochrome c", "ko:K00401".

Extract all genes that have a gene qualifier that equals 'omcX'

::

    biosqlx export sequence -q gene -v omcX -o feat-prot

Extract all genes that have the EC number 2.1.1.1 but only if they are
from Archaea

::

    biosqlx export sequence -o feat-prot -t Archaea -q EC_number -v 2.1.1.1

If you want to extract sequences using a KEGG ortholog number then you
need to use ``db_xref`` as the qualifier and have ``ko:`` before the
KEGG accession for the particular gene. This is an unintuitive syntax
but is required as KEGG ortholog information is stored specially in the
database.

::

    biosqlx export sequence -q db_xref -v ko:K00399 -o feat-prot

splitting output into separate files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default all of the output will be printed to stdout, which can then
be redirected to a file. However it's also possible to output
information to files for individual species using the
``--split-species`` option.

the following will create individual fasta files for all of the species
that belong to desulfovibrio in the database

::

    biosqlx export sequence -t Desulfovibrio --split-species

``biosqlx export taxonomy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This subcommand will allow you to see the organsims currently in the
database. The data is presented as either a hierarchical tree or as 
a semicolon separated list of taxonomy strings. The ``--root`` option
allows you to output only those organisms that fall under the named
taxonomy.

``biosqlx export namespace``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will print the information about the various samples/namespaces
present in the database.

``biosqlx modify``
------------------


``biosqlx modify annotation``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subcommand will add or modify an annotation to a seqfeature (gene).
You provide a **tab separated** input file that describes the
annotations to add, where the first row *must* be a header that
describes the qualifiers to add and one of the columns *must* uniquely
identify a seqfeature. The name of this column must be given using the
``--key`` agrument on the command line

Lets look at an example of what "qualifiers" are and how they could be
represented in the input file. Below is a excerpt from a genbank file
that shows all of the annotations for a particular protein. The
qualifiers of the gene are shown on the lines that begin with a "/"
character and come before the "=" character (eg. gene, EC\_number).

::

    CDS             complement(6523..7818)
                    /gene="proS"
                    /locus_tag="KQ51_00006"
                    /EC_number="6.1.1.15"
                    /inference="ab initio prediction:Prodigal:2.60"
                    /inference="similar to AA sequence:UniProtKB:A6U7Z3"
                    /codon_start=1
                    /transl_table=11
                    /product="Proline--tRNA ligase"
                    /protein_id="AIO17910.1"
                    /db_xref="GI:685629398"

This could be mapped onto a row of the input file as follows

::

    locus_tag   gene    EC_number   product
    KQ51_00006  proS    6.1.1.15    Proline--tRNA ligase

In this case the key column is ``locus_tag`` with the value
``KQ51_00006`` and then add the values to that gene for the given
qualifier. The key column, whatever tag it is, must be unique amongst
all genes in the database. Good qualifiers to use would be
``locus_tag``, ``ID`` or ``protein_id`` as they are often unique.
**However**, none of these qualifiers are *guaranteed* to be unique in
our database, so be careful. The database itself has an ID called a
``seqfeature_id`` that *is* guaranteed to be unique, so if you know the
seqfeature\_id then use that. If the key column is the seqfeature\_id,
then you must provide the ``-s`` option to the script.

When adding annotations to a gene, the default behavior is to add a
second annotation to a gene if one already exists for that qualifier.
This may not be what you want, say if the original annotation is
incorrect; in this case use the ``--replace`` flag to the script.

::

    biosqlx modify annotation -i annotations.tsv --key seqfeature_id
    biosqlx modify annotation -i annotations.tsv --replace --key locus_tag


``biosqlx modify taxonomy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Add, change, or remove taxonomy IDs for sequences or the taxonomy
tree itself.

Moving a taxon underneath a new parent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``--move`` option allows you to change the structure of the tree by
letting you move a taxon underneath a new parent node. Two arguments
need to be given, first the child taxon (the one you want to move)
and then the parent taxon::

    biosqlx modify taxonomy --move "Methanosaeta harundinacea" Methanothrix

The above command moves the species *Methanosaeta harundinacea* underneath
the genus *Methanothrix*. Both of these taxonomy names must exist already
in the database for the opperation to take place. Notice the quotes
surronding *Methanosaeta harundinacea*, they are required whenever a
taxon name is more than one space separated word.

It is also possible to use an NCBI taxonomy ID instead of a taxon name
for either the child or parent taxons. The example above could be written
in any of the following ways::

    biosqlx modify taxonomy --move 2223 2222
    biosqlx modify taxonomy --move "Methanosaeta harundinacea" 2222
    biosqlx modify taxonomy --move 2223 Methanothrix

This requires that these taxons have the NCBI taxon ID associated with them.

``biosqlx add``
---------------


``biosqlx add sequence``
~~~~~~~~~~~~~~~~~~~~~~~~

This is the main way to add in new datasets to the database. You'll
need to have the sequences at least run through an ORF caller, such as
Prodigal, to add them into the database. Sequences can either be given
as a genbank formatted file, using the ``-G`` option or be provided as a
fasta plus gff files, using the ``-f`` and ``-g`` options. At this time
plain fasta files without any ORFs called are not supported.

Specifiying an NCBI taxonomy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Genbank files can contain information about the taxonomy of the organism,
which can be used to populate the taxonomy in the database. By default,
minimal taxonomy information is added to the database when a new sequence
is added from a taxon that is not currently in the database. For example,
if you were to load the genome of *Methanosaeta harundinacea* into the
database, when no other archaea were present, then only the organism name
will be stored. However, by specifying ``-t``, the taxonomy information
will be downloaded from NCBI and the full taxonomy tree will be populated
for the organism.

The second option is to combine ``-t`` with ``-T`` to provide a NCBI
taxon ID on the commandline. This is useful if the input file doesn't
contain taxonomic information for the organism. For example, specifying
a gff file with ``-g`` that does not contain taxonomy information. Or
alternatively it can be used to overwrite the taxonomy information given
in the input files.

Specifying new taxonomies
^^^^^^^^^^^^^^^^^^^^^^^^^
Sometimes you may need to add novel organisms that are not currently in
the NCBI taxonomy database. In this situation you can specify new taxons
on the commandline using the format ``<taxon_name>:<taxon_rank>``, where
``<taxon_name>`` is the new name that you with to add and ``<taxon_rank>``
is a recognized taxonomic rank, such as "kingdom", "phylum", "genus",
"species". It is also possible to specify multiple taxons in order of
increasing specificity, for example::

    biosqlx add sequence -T 2 -t -G GCA_000830255.1.gb Epsilonbacterota:phylum Campylobacteria:class Campylobacterales:order Thiovulaceae:family PC08-66:genus "Sulfuricurvum sp. PC08-66:species"

Notice above that the new taxonomy is listed in increasing specificity
(phylum, class, order, family, genus, species), and the quotes around
the species name, since the name contains space characters. The ``-T 2``
in this example means that the novel taxons listed on the commandline
begin (or are children) of the NCBI taxon ID (in this case 2 equals
bacteria). Any NCBI taxon can be given and the new taxons will be
children::

    biosqlx add sequence -T 94695 -t -G ANME_genome.gb ANME-2ab:family ANME-2b:genus "ANME sp. NewGenome:species"

In this example the NCBI taxon, 94695 is the order *Methanosarcinales* and
the new taxons specified on the commandline give a novel family, genus and
species. The new taxons given on the commandline are also checked against
the database when adding new sequences, so the following will also work::

    biosqlx add sequence -T 94695 -t -G ANME_genome2.gb ANME-2ab:family ANME-2a:genus "ANME sp. AnotherNewGenome:species"

In this case the novel family, ANME-2ab, is already in the database,
from the previous example, and so is not added again. The genus and
species are novel and will be added into the database.

Be careful about how you arrange the new taxons on the commandline;
they must be in the correct order as no checking is performed on the
``<taxon_rank>`` itself. It's possible to specify something like
the following: ``"ANME sp. AnotherNewGenome:species" ANME-2a:genus
ANME-2ab:family`` which will not cause an error and actually produce
the following tree::

    #This is the opposite of what was intended 
    ├── ANME sp. AnotherNewGenome (species)
    │   ├── ANME-2a (genus)
    │   │   └── ANME-2ab (family)


