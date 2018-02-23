# Introduction
This is a repository of sql tables and views and scripts that could be added onto an
existing BioSQL database. There are a number of scripts that should hopefully work on
any BioSQL database and have been developed for easy access to various information.


In addition there tables that were made solely for sqlite3
and are primarily related to adding functionality for performing proteomics
experiments. I have not tried to write the proteomics tables in a porable way, they were designed to model the data created by the SIPROS program.


## Python API    

### `extract_sequences_from_taxonomy.py`
This script will extract from the database all of the sequences associated with
a particular taxonomy. The input is either an NCBI taxonomy ID or the complete
taxonomic name. The output, specified using the `-o` option can be one of fasta,
gb (genbank), feat-nucl (genes output as fasta nucleotide) or feat-prot (genes
output as fasta amino acids)

#### Examples
Extract translated gene sequences from all of the Deltaproteobacteria in the database
```
extract_sequences_from_taxonomy.py -u orphanlab -o feat-prot -d biosqldb Deltaproteobacteria
```

Extract all of the ANME-2a sequences as a genbank file specifying the host for the
database. This option is required when running the script on a computer other than
ocean (like your own laptop)
```
extract_sequences_from_taxonomy.py -u orphanlab -o gb -H ocean.gps.caltech.edu -d biosqldb ANME-2a
```

Extract all of the sequences from the database that match the NCBI taxonomy ID, 872.
In this case that refers to the [Desulfovibrio](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=872&lvl=3&lin=f&keep=1&srchmode=1&unlock).
```
extract_sequences_from_taxonomy.py -u orphanlab biosqldb 872
```

When using the taxonomic name has the limitation that it must be a complete name
so for example if you wanted sequences from a species then you must give both
the genus and species name. In the below example, if you just type 'vulgaris',
no results will be returned.
```
extract_sequences_from_taxonomy.py -u orphanlab -o feat-prot biosqldb 'Desulfovibrio vulgaris'
```
Note the single quotes around the organism name! This is important so that the
name is properly input into the script. Any names that are multiple words need
to have the quotes surrounding them.

### `add_annotation_to_protein.py`
This script will add an annotation to a seqfeature (gene). You provide a **tab
separated** input file that describes the annotations to add, where the first
row *must* be a header that describes the qualifiers to add and one of the columns
*must* uniquely identify a seqfeature. The name of this column must be given using
the `--key` agrument on the command line

Lets look at an example of what "qualifiers" are and how they could be represented
in the input file. Below is a excerpt from a genbank file that shows all of the
annotations for a particular protein. The qualifiers of the gene are shown on
the lines that begin with a "/" character and come before the "=" character
(eg. gene, EC_number).

```
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
```

This could be mapped onto a row of the input file as follows
```
locus_tag   gene    EC_number   product
KQ51_00006  proS    6.1.1.15    Proline--tRNA ligase
```
In this case the key column is `locus_tag` with the value `KQ51_00006`
and then add the values to that gene for the given qualifier. The key column,
whatever tag it is, must be unique amongst all genes in the database. Good qualifiers
to use would be `locus_tag`, `ID` or `protein_id` as they are often unique. **However**,
none of these qualifiers are *guaranteed* to be unique in our database, so be careful.
The database itself has an ID called a `seqfeature_id` that *is* guaranteed to be
unique, so if you know the seqfeature_id then use that. If the key column is
the seqfeature_id, then you must provide the `-s` option to the script.

When adding annotations to a gene, the default behavior is to add a second
annotation to a gene if one already exists for that qualifier. This may not be
what you want, say if the original annotation is incorrect; in this case use the
`--replace` flag to the script.

#### Examples

 ```
 add_annotation_to_protein.py -u orphanlab -i annotations.tsv -s -d biosqldb --key seqfeature_id
 ```

 ```
 add_annotation_to_protein.py -u orphanlab -i annotations.tsv --replace -d biosqldb --key locus_tag
 ```

### `extract_sequences_using_qv.py`
This script will extract all of the genes that meet certain criteria based on
the values of their qualifiers. This is a strict evaluation of equality so everything
must be spelled correctly etc. Unfortunately this also means that you cannot do
other operations like asking for all features that are greater than or less than
some value.

#### Examples

Extract all genes that have a gene qualifier that equals 'omcX'
```
extract_sequences_using_qv.py -u orphanlab -d biosqldb gene omcX
```

Extract all genes that have the EC number 2.1.1.1 but only if they are from the
5133 megahit metagenome
```
extract_sequences_using_qv.py -u orphanlab -d biosqldb -D 5133_megahit EC_number 2.1.1.1
```

If you want to extract sequences using a KEGG ortholog number then you need to
use `db_xref` as the qualifier and have `ko:` before the KEGG accession for the
particular gene. This is an unintuitive syntax but is required as KEGG ortholog
information is stored specially in the database.  
```
extract_sequences_using_qv.py -u orphanlab -d biosqldb db_xref ko:K00399
```

### `extract_seqfeatures.py`
This is a lower level script that will extract genes using the seqfeature_id, which
is the database specific identifier that is guaranteed to be unique. Supply an
input file containing seqfeature_ids. Will print fasta files of the genes.

### `get_annotations_for_seqfeature.py`
This script prints all of the qualifier value information for each gene as a
csv file. You must provide an input file containing seqfeature_ids

### `dump_biodatabae.py`
output all of the genes in a particular biodatabase   
