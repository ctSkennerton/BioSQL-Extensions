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
------------

There is a main script ``biosqlx`` that serves as the main extry point
for quickly accessing information in the database. There are a number of
subcommands that allow you to add, modify or export information from the
database. For example ``biosqlx add sequence`` allows you to add new
sequence data into the database, while ``biosqlx export sequence`` lets
you extract data out.

connecting to the database
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
^^^^^^^^^^^

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

``biosqlx export sequence``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subcommand allows you to slice and dice the sequence data using
multiple filtering metrics, including the assigned taxonomy or the
sequence feature annotations.

Output formats
^^^^^^^^^^^^^^

There are currently five supported output formats: genbank, fasta,
feat-prot, feat-nucl, csv. Genbank outputs a genbank formatted file for
the sequence records, these are likely contigs or whole chromosomes.
Fasta outputs the sequence information in fasta format. feat-prot and
feat-nucl output the sequence features, like protein sequences, tRNAs,
rRNAs. feat-prot outputs the translated features, so by default will
only output the proteins. feat-nucl will output untranslated features
and by default outputs the proteins, tRNAs and rRNAs. The csv format
will output the annotations for all of the features.

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
