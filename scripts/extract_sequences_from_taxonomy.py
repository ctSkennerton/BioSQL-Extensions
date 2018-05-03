#!/usr/bin/env python
from __future__ import print_function
import sys
from getpass import getpass, getuser
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeqRecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from common import generate_placeholders, chunks, extract_feature_sql, standard_options, get_seqfeature_ids_for_bioseqs
from get_annotations_for_seqfeature import print_feature_qv_csv

def extract_feature(dbrec, output_format, fp, wanted_types=['CDS','rRNA', 'tRNA'], id_tag=None):

    for feature in dbrec.features:
        # only print the proteins
        if feature.type not in wanted_types:
            continue

        if 'pseudo' in feature.qualifiers:
            continue

        feat_extract = feature.extract(dbrec.seq.toseq())
        codon_start = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
        feat_extract = feat_extract[codon_start:]
        if output_format == 'feat-prot':
            feat_extract = feat_extract.translate()

        if id_tag:
            try:
                seqid = feature.qualifiers[id_tag][0]
            except KeyError:
                print("WARNING: cannot find {} tag for seqfeature {}".format(id_tag, feature._seqfeature_id), file=sys.stderr)
        else:
            seqid = feature._seqfeature_id

        description = ''
        try:
            description += ' '.join(feature.qualifiers['product'])
        except KeyError:
            pass
        else:
            description += ' '

        try:
            description += ' '.join(feature.qualifiers['gene'])
        except KeyError:
            pass

        feat_extract = SeqRecord(feat_extract, id=str(seqid), description=description)
        SeqIO.write(feat_extract, fp, 'fasta')

def main(args):
    server = BioSeqDatabase.open_database(driver=args.driver, db=args.database, user=args.user, host=args.host, passwd=args.password)

    tax_name = False
    try:
        ncbi_tax = int(args.taxid)
    except ValueError:
        tax_name = True

    if not tax_name:
        print("interpreting as an NCBI taxon ID...", file=sys.stderr)
        taxon_id_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.ncbi_taxon_id  = %s AND include.right_value = include.left_value + 1)"

        rows = server.adaptor.execute_and_fetchall(taxon_id_lookup_sql, (ncbi_tax,))
    else:
        print("interpreting as a taxon name...", file=sys.stderr)
        taxon_name_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name "\
                "WHERE name like %s) AND include.right_value = include.left_value + 1)"
        rows = server.adaptor.execute_and_fetchall(taxon_name_lookup_sql, (args.taxid,))

    if args.feature_type is not None:
        types = args.feature_type
    elif args.output_format == 'feat-prot':
        types = ['CDS']
    elif args.output_format == 'feat-nucl':
        types = ['CDS', 'rRNA', 'tRNA']

    if len(rows) == 0:
        print("\nThere does not appear to be any sequences associated with\n"
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
    if args.split_species:
        taxon_file_mapping = {}
        for k, v in dbids.items():
            tname = server.adaptor.execute_and_fetch_col0("SELECT name from taxon_name where taxon_id = %s and name_class = %s", (v,'scientific name'))[0]
            tname = tname.replace(' ', '_')
            if args.output_format == 'gb':
                tname += '.gb'
            elif args.output_format == 'feat-prot':
                tname += '.faa'
            elif args.output_format == 'csv':
                tname += '.csv'
            else:
                tname += '.fna'
            files[v] = tname
            taxid_to_dbids.setdefault(v, []).append(k)


    if args.split_species:
        # got to save all of the records before printing them out
        outdata = {}
        for taxid, dbid_list in taxid_to_dbids.items():
            if args.output_format == 'csv':
                with open(files[taxid], 'w') as fp:
                    print_feature_qv_csv(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbid_list]), fp)
            else:
                for dbid, dbname in dbid_list:
                    db = server[dbname]
                    seq_rec = db[dbid]
                    outdata.setdefault(taxid, []).append(seq_rec)

        for taxid, dbrecs in outdata.items():
            with open(files[taxid], 'w') as fp:
                if 'feat' in args.output_format:
                    for dbrec in dbrecs:
                        extract_feature(dbrec, args.output_format, fp)
                elif 'csv' != args.output_format:
                    SeqIO.write(dbrecs, fp, args.output_format)

    else:
        if args.output_format == 'feat-prot':
            extract_feature_sql(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]),type=types, translate=True )
        elif args.output_format == 'feat-nucl':
            extract_feature_sql(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]), type=types)
        elif args.output_format == 'csv':
            print_feature_qv_csv(server, get_seqfeature_ids_for_bioseqs(server, [x[0] for x in dbids.keys()]))
        else:
            for (dbid, dbname), taxid in dbids.items():
                db = server[dbname]
                try:
                    dbrec = db[dbid]
                    SeqIO.write(dbrec, sys.stdout, args.output_format)
                except KeyError:
                    pass


if __name__ == "__main__":
    parser = standard_options(description="This script will extract from the database all of the sequences associated with a particular taxonomy. The input is either an NCBI taxonomy ID or the complete taxonomic name.")
    parser.add_argument('-o', '--output_format', help='output format of the selected sequences. Choices: fasta - fasta file of the contigs; gb - genbank file of the sequences; feat-prot - fasta file containing the translated coding sequences; feat-nucl - fasta file containing the untranslated coding sequences, tRNAs and rRNAs; csv - csv file of annotations for the features', choices=['fasta', 'gb', 'feat-prot', 'feat-nucl', 'csv'], default='fasta')
    parser.add_argument('taxid', help='supply a ncbi taxonomy id that will be extracted. If an integer is supplied it will be interpreted as an NCBI taxonomy id; otherwise it will be interpreted as part of a taxonomy name (e.g. Proteobacteria)', default=None)
    parser.add_argument('-s', '--split_species', help='when there are multiple species to be returned, split them into separate files, based on their name, instead of printing to stdout', default=False, action='store_true')
    parser.add_argument('-t', '--feature-type', help='restrict the results to feature type e.g. rRNA, tRNA, CDS. This option can be specified multiple times for multiple types', default=None, action='append')
    args = parser.parse_args()
    if args.password is None:
        args.password = getpass("Please enter the password for user " + \
                args.user + " on database " + args.database)
    main(args)

