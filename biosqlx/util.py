import sys
import csv
from Bio.Seq import reverse_complement, translate as bio_translate
class OrthologData(object):
    """
    """
    def __init__(self, term_id, term_name):
        self.term_id = term_id
        self.term_name = term_name
        self.product = ""
        self.ec = []
        self.dbxref_id = None
        self.gene = []
        self.genomes = {}

class DbError(Exception):
    pass

class DbInputError(DbError):
    def __init__(self, message):
        self.message = message

def generate_placeholders(l):
    placeholder= ['%s'] # use ? For SQLite. See DBAPI paramstyle.
    return ', '.join(placeholder * l)


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def get_organism_name(server, bioid):
    sql = "SELECT name from taxon_name where taxon_id = \
            (SELECT taxon_id from bioentry where bioentry_id = %s) \
            and name_class = 'scientific name'"
    return server.adaptor.execute_and_fetchall(sql, (bioid,))[0][0]

def standard_options(**kwargs):
    import argparse
    from getpass import getuser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, **kwargs)
    parser.add_argument('-d', '--database', help='name of premade biosql database', required=True)
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host',
            default=5432)
    parser.add_argument('-u', '--user', help='database user name',
            default=getuser())
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to',
            default='localhost')

    return parser


def get_bioentries_from_taxonomy(server, taxid):
    tax_name = False
    try:
        ncbi_tax = int(taxid)
    except ValueError:
        tax_name = True

    if not tax_name:
        #print("interpreting as an NCBI taxon ID...", file=sys.stderr)
        taxon_id_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.ncbi_taxon_id  = %s AND include.right_value = include.left_value + 1)"

        rows = server.adaptor.execute_and_fetchall(taxon_id_lookup_sql, (ncbi_tax,))
    else:
        #print("interpreting as a taxon name...", file=sys.stderr)
        taxon_name_lookup_sql = "SELECT bioentry_id, taxon_id, biodatabase.name FROM bioentry JOIN "\
                "biodatabase USING(biodatabase_id) WHERE taxon_id IN "\
                "(SELECT DISTINCT include.taxon_id FROM taxon "\
                "INNER JOIN taxon as include ON (include.left_value "\
                "BETWEEN taxon.left_value AND taxon.right_value) "\
                "WHERE taxon.taxon_id IN (SELECT taxon_id FROM taxon_name "\
                "WHERE name like %s) AND include.right_value = include.left_value + 1)"
        rows = server.adaptor.execute_and_fetchall(taxon_name_lookup_sql, (taxid,))

    return rows


def print_feature_qv_csv(server, sfids, outfile=sys.stdout):
    """raw sql extraction of fasta seqfeatures
    """
    def dbxref_dict(server, seqfeature_ids):
        db_qv = {}
        for feat_chunk in chunks(seqfeature_ids, 900):
            #sql = "SELECT s.seqfeature_id, d.dbname || ':' || d.accession AS kegg_id "\
            #        "FROM seqfeature_dbxref s "\
            #        "JOIN dbxref d USING(dbxref_id) "\
            #        "WHERE s.seqfeature_id IN ({})".format(generate_placeholders(len(feat_chunk)))
            sql = "SELECT s.seqfeature_id, d.dbname, d.accession, t.name, dqv.value "\
                    "FROM seqfeature_dbxref s "\
                    "JOIN dbxref d USING(dbxref_id) "\
                    "LEFT JOIN dbxref_qualifier_value dqv USING(dbxref_id) "\
                    "LEFT JOIN term t USING(term_id) "\
                    "WHERE s.seqfeature_id IN ({}) "\
                    "ORDER BY s.seqfeature_id, d.dbname, s.rank".format(generate_placeholders(len(feat_chunk)))
            for seqfeature_id, dbname, dbxref, name, value in server.adaptor.execute_and_fetchall(sql, tuple(feat_chunk)):
            #for seqfeature_id, dbxref in server.adaptor.execute_and_fetchall(sql, tuple(feat_chunk)):
                try:
                    db_qv[seqfeature_id][dbname] = dbxref
                except KeyError:
                    db_qv[seqfeature_id] = {}
                    db_qv[seqfeature_id][dbname] = dbxref

                if name:
                    db_qv[seqfeature_id][name] = value
        return db_qv


    def qv_dict(server, seqfeature_ids):
        qv = {}
        for feat_chunk in chunks(seqfeature_ids, 900):
            feat_chunk2 = tuple(feat_chunk)
            qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv JOIN term t ON t.term_id = qv.term_id WHERE seqfeature_id IN ({})'.format(generate_placeholders(len(feat_chunk)))

            taxonomy_sql = 'SELECT seqfeature_id, bioentry.name, biodatabase.name, lineage.lineage FROM seqfeature JOIN bioentry USING(bioentry_id) JOIN biodatabase USING(biodatabase_id) LEFT JOIN lineage ON taxon_id = lineage.id WHERE seqfeature_id IN ({})'.format(generate_placeholders(len(feat_chunk)))
            for seqfeature_id, contig, namespace, lineage in server.adaptor.execute_and_fetchall(taxonomy_sql, feat_chunk2):
                try:
                    if lineage:
                        qv[seqfeature_id]['taxonomy'] = lineage
                        qv[seqfeature_id]['organism'] = lineage.split(';')[-1]

                    qv[seqfeature_id]['bioentry'] = contig
                    qv[seqfeature_id]['sample'] = namespace
                except KeyError:
                    qv[seqfeature_id] = {}
                    if lineage:
                        qv[seqfeature_id]['taxonomy'] = lineage
                        qv[seqfeature_id]['organism'] = lineage.split(';')[-1]

                    qv[seqfeature_id]['bioentry'] = contig
                    qv[seqfeature_id]['sample'] = namespace

            for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, feat_chunk2):
                if not name:
                    continue
                try:
                    qv[seqfeature_id][name] = value
                except KeyError:
                    qv[seqfeature_id] = {}
                    qv[seqfeature_id][name] = value

        return qv


    qv = qv_dict(server, sfids)
    dbxr = dbxref_dict(server, sfids)
    for sf, data in dbxr.items():
        for n, v in data.items():
            qv[sf][n] = v

    columns = set()
    for sf, data in qv.items():
        columns |= set(data.keys())

    columns = sorted(columns)
    writer = csv.writer(outfile)
    writer.writerow(['seqfeature'] + columns)
    for sf, data in qv.items():
        row = [sf]
        for i in columns:
            row.append(data.get(i, None))
        writer.writerow(row)

def get_seqfeature_id_from_qv(db, qualifier, value, biodatabase_id=None):
    rows = db.adaptor.execute_and_fetchall('select term_id from term join ontology using(ontology_id) where term.name = %s and ontology.name = \'Annotation Tags\'', (qualifier,))
    if len(rows) != 1:
        raise ValueError("The qualifier is not unique")

    term_id = rows[0][0]
    sql = r'select seqfeature_id from seqfeature_qualifier_value join seqfeature using(seqfeature_id) join bioentry using(bioentry_id) where seqfeature_qualifier_value.term_id = %s and value = %s'
    if biodatabase_id is not None:
        sql += ' and biodatabase_id = %s'

    if biodatabase_id is not None:
        rows = db.adaptor.execute_and_fetchall(sql, (term_id, value, db.dbid))
    else:
        rows = db.adaptor.execute_and_fetchall(sql, (term_id, value))

    if len(rows) > 1:
        raise ValueError("There is more than one seqfeature associated with qualifier={} and value={}".format(qualifier, value))
    elif len(rows) == 0:
        raise ValueError("There are no seqfeature associated with qualifier={} and value={}".format(qualifier, value))

    return rows[0][0]


def get_seqfeature_ids_from_qv(db, qualifier, value, namespace=None, fuzzy=False):
    ''' Return a set of seqfeatures based on a qualifier, value and database.

    :param db: Ther server instance
    :param qualifier: text of the qualifier name to search for. eg. product, gene
    :param value: text of the value of the qualifier. eg. Hydrogenase, mcrA
    :param namespace: The name of the sub-database to restrict matches to.
        Corresponds to biodatabase.name in the SQL schema.
    :param fuzzy: Set to True to use the LIKE command matching rather than exact matching.
    :returns: list of seqfeature_ids (integers)
    :raises DbInputError
    '''
    if qualifier == 'db_xref':
        # need to handle the special instance of qualifiers that refer to other databases
        try:
            dbname, accession = value.split(':')
            if namespace is not None:
                sql = 'SELECT seqfeature_id FROM seqfeature_dbxref ' \
                      'JOIN dbxref ON dbxref.dbxref_id = seqfeature_dbxref.dbxref_id ' \
                      'JOIN seqfeature ON seqfeature_dbxref.seqfeature_id = seqfeature.seqfeature_id ' \
                      'JOIN bioentry ON bioentry.bioentry_id = seqfeature.bioentry_id ' \
                      'JOIN biodatabase ON biodatabase.biodatabase_id = bioentry.biodatabase_id ' \
                      'WHERE dbxref.dbname = %s AND dbxref.accession = %s AND biodatabase.name = %s'
                col0 = db.adaptor.execute_and_fetch_col0(sql, (dbname, accession, namespace))
            else:
                sql = 'SELECT seqfeature_id FROM seqfeature_dbxref ' \
                      'JOIN dbxref ON dbxref.dbxref_id = seqfeature_dbxref.dbxref_id ' \
                      'WHERE dbxref.dbname = %s AND dbxref.accession = %s'
                col0 = db.adaptor.execute_and_fetch_col0(sql, (dbname, accession))
        except ValueError:
            raise DbInputError('''Error: value does not contain both a database name and value
Hint: When using db_xref as the input qualifier, the value must contain two terms separated
by a colon (:) character, for example ko:K03388, where the first part is the database name
and the second part is the value in that crossreferenced database. The offending value is: ''' + value)


    else:
        if fuzzy:
            if namespace is not None:
                sql = 'SELECT qv.seqfeature_id FROM seqfeature_qualifier_value qv ' \
                      'JOIN seqfeature s ON qv.seqfeature_id=s.seqfeature_id ' \
                      'JOIN bioentry b ON b.bioentry_id=s.bioentry_id ' \
                      'JOIN term t ON t.term_id=qv.term_id ' \
                      'JOIN biodatabase d ON d.biodatabase_id=b.biodatabase_id ' \
                      'WHERE t.name = %s AND qv.value LIKE %s AND d.name = %s'
            else:
                sql = 'SELECT seqfeature_id FROM seqfeature_qualifier_value ' \
                      'JOIN term USING(term_id) WHERE term.name = %s AND value LIKE %s'
        else:
            if namespace is not None:
                sql = 'SELECT qv.seqfeature_id FROM seqfeature_qualifier_value qv ' \
                      'JOIN seqfeature s ON qv.seqfeature_id=s.seqfeature_id ' \
                      'JOIN bioentry b ON b.bioentry_id=s.bioentry_id ' \
                      'JOIN term t ON t.term_id=qv.term_id ' \
                      'JOIN biodatabase d ON d.biodatabase_id=b.biodatabase_id ' \
                      'WHERE t.name = %s AND qv.value = %s AND d.name = %s'
            else:
                sql = 'SELECT seqfeature_id FROM seqfeature_qualifier_value ' \
                      'JOIN term USING(term_id) WHERE term.name = %s AND value = %s'


        if namespace is not None:
            col0 = db.adaptor.execute_and_fetch_col0(sql, (qualifier, value, biodatabase_id))
        else:
            col0 = db.adaptor.execute_and_fetch_col0(sql, (qualifier, value))

    if len(col0) == 0:
        raise ValueError("There are no seqfeature associated with "
                         "qualifier={} and value={}".format(qualifier, value))

    return col0

def get_seqfeature_from_input(server, input_ids, type='ID'):
    ret = []
    for c in chunks(input_ids, 900):
        sql = "SELECT seqfeature_id FROM seqfeature_qualifier_value where term_id = \
                (SELECT term_id FROM term WHERE name = %s) AND value IN ({})"
        sql = sql.format(generate_placeholders(len(c)))
        for row in server.adaptor.execute_and_fetchall(sql, tuple([type] + c)):
            ret.append(row[0])

    return ret

def get_bioseqid_for_seqfeature(server, ids):
    '''
    '''
    bioentry_ids = []
    for c in chunks(ids, 900):
        sql = "SELECT d.name, s.bioentry_id, s.seqfeature_id, b.taxon_id FROM seqfeature s  \
                JOIN bioentry b ON s.bioentry_id = b.bioentry_id JOIN biodatabase d ON b.biodatabase_id = d.biodatabase_id WHERE \
                s.seqfeature_id IN ({})".format(generate_placeholders(len(c)))
        for row in server.adaptor.execute_and_fetchall(sql, tuple(c)):
            bioentry_ids.append(row)

    return bioentry_ids

def get_seqfeature_ids_for_bioseqs(server, ids):
    seqfeature_ids = []
    for c in chunks(ids, 900):
        sql = "SELECT seqfeature_id FROM seqfeature WHERE bioentry_id IN ({})".format(generate_placeholders(len(c)))
        for row in server.adaptor.execute_and_fetchall(sql, tuple(c)):
            seqfeature_ids.append(row[0])
    return seqfeature_ids

def extract_feature_sql(server, seqfeature_ids, type=['CDS', 'rRNA', 'tRNA'], qualifier=['ID','locus_tag'], translate=False, file=sys.stdout):
    """raw sql extraction of fasta seqfeatures
    """
    for chunk in chunks(seqfeature_ids, 900):
        sql = "SELECT f.seqfeature_id AS gid, \
                      fl.strand,\
                      substr(s.seq, COALESCE(fl.start_pos, 1), (COALESCE(fl.end_pos, s.length) - COALESCE(fl.start_pos, 1))+1) AS subseq\
               FROM   seqfeature f \
               JOIN   term t ON f.type_term_id=t.term_id \
               JOIN   location fl USING(seqfeature_id) \
               JOIN   biosequence s USING(bioentry_id)\
               WHERE  t.name IN ({}) AND f.seqfeature_id IN ({})".format(generate_placeholders(len(type)), generate_placeholders(len(chunk)))

        features = server.adaptor.execute_and_fetchall(sql, tuple(type + chunk) )
        results = {}
        for sfid, strand, seq in features:
            if seq is None:
                # at the moment partial features are stored with a null at either the start_pos or end_pos
                # in the location table in the database. There isn't a good way at the moment
                print("seqfeature {} could not be extracted, problem with the location".format(sfid), file=sys.stderr)
            else:
                results[sfid] = (strand, seq)

        qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv, term t WHERE seqfeature_id IN ({}) AND t.term_id = qv.term_id'.format(generate_placeholders(len(chunk)))
        qv = {}
        for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, tuple(chunk)):
            try:
                qv[seqfeature_id][name] = value
            except KeyError:
                qv[seqfeature_id] = {}
                qv[seqfeature_id][name] = value

        tax_name_select_sql = 'SELECT seqfeature.seqfeature_id, taxon_name.name FROM seqfeature JOIN bioentry ON seqfeature.bioentry_id = bioentry.bioentry_id JOIN taxon_name ON bioentry.taxon_id = taxon_name.taxon_id WHERE seqfeature.seqfeature_id IN ({}) AND taxon_name.name_class = \'scientific name\''.format(generate_placeholders(len(chunk)))
        tax = {}
        for seqfeature_id, name in server.adaptor.execute_and_fetchall(tax_name_select_sql, tuple(chunk)):
            tax[seqfeature_id] = name

        for seqfeature_id, (strand, seq) in results.items():
            # remove any pseudo genes
            if 'pseudo' in qv[seqfeature_id]:
                continue

            name = str(seqfeature_id)
            for q in qualifier:
                try:
                    name += ' ' + qv[seqfeature_id][q]

                except KeyError:
                    pass

            if strand == -1:
                try:
                    seq = reverse_complement(results[seqfeature_id][1])
                except TypeError as e:
                    raise TypeError("failed to retieve sequence for {}".format(seqfeature_id))

            try:
                codon_start = int(qv[seqfeature_id]['codon_start'][0]) - 1
            except KeyError:
                try:
                    codon_start = int(qv[seqfeature_id]['phase'][0])
                except KeyError:
                    codon_start = 0

            seq = seq[codon_start:]
            if translate:
                seq = bio_translate(seq)
            try:
                name += ' ' + qv[seqfeature_id]['product']
            except KeyError:
                pass

            try:
                name += ' [' + tax[seqfeature_id] + ']'
            except KeyError:
                pass

            print(">{}\n{}".format(name, seq), file=file)

def get_kegg_id_from_name(server, orthology, brite='KEGG'):
    ''' get the ID for the wanted orthology
        server      the database connection object
        brite       the particular type of orthology like KEGG or KEGG_modules
        orthology   the actual term to search for

        returns the internal database ID for the ortholog
    '''

    sql = 'select term_id from term where name like ?  and ontology_id = (select ontology_id from ontology where name = %s)'
    rows = server.adaptor.execute_and_fetchall(sql, ('%'+orthology+'%', brite))
    if len(rows) != 1:
        raise DbError('That orthology either doesn\'t exist or isn\'t unique')

    return rows[0][0]


def get_kegg_data_from_id(server, orthology_id, brite="KEGG"):
    ''' Get the dbxref data for an ortholgy
        server          the database server object
        orthology_id    the term_id for the orthology
        brite           The kegg ontology to use -> KEGG or KEGG_modules

        returns a dict of dbxref_id -> OrthologyData objects
    '''
    # get a list of all orthologs underneath the wanted module
    sql = 'select object_term_id, term.name, term.identifier, distance from term_path join term on object_term_id = term.term_id where subject_term_id = %s and term_path.ontology_id = (select ontology_id from ontology where name = %s)'
    rows = server.adaptor.execute_and_fetchall(sql, (orthology_id, brite))

    term_names = {}
    for x in rows:
        if re.match(r'K\d+', x[1]):
            term_names[x[0]] = x[1]

    #print("found {} orthologs".format(len(term_names)), file=sys.stderr )

    term_ids = term_names.keys()  #[x[0] for x in rows if x[3] == max_distance]
    dbxref_ids = []
    kegg_orthology_data = {}

    for chunk in chunks(term_ids, 990):
        sql = 'select dbxref_id, term_id from term_dbxref where term_id in (%s)'
        sql = sql % generate_placeholders(len(chunk))
        rows = server.adaptor.execute_and_fetchall(sql, tuple(chunk))
        for row in rows:
            kegg_orthology_data[row[0]] = OrthologData(row[1], term_names[row[1]])
            kegg_orthology_data[row[0]].dbxref_id = row[0]
            dbxref_ids.append(row[0])


    # Now get some of the metadata for the dbxref like the name
    sql = 'select * from dbxref_qualifier_value where term_id = (select term_id from term where name = \'product\' ) and dbxref_id in (%s)'
    for chunk in chunks(dbxref_ids, 990):
        sql = sql % generate_placeholders(len(chunk))
        rows = server.adaptor.execute_and_fetchall(sql, tuple(chunk))
        for row in rows:
            kegg_orthology_data[row[0]].product = row[3]

    return kegg_orthology_data


def filter_seqfeature_ids_by_taxonomy(server, seqfeature_ids, taxid):
    wanted_bioentries = get_bioentries_from_taxonomy(server, taxid)
    seqfeature_bioentries = get_bioseqid_for_seqfeature(server, seqfeature_ids)
    seqfeature_bioentry_map = {}
    for row in seqfeature_bioentries:
        try:
            seqfeature_bioentry_map[row[1]].append(row[2])
        except KeyError:
            seqfeature_bioentry_map[row[1]] = [row[2]]

    final_seqfeatures = []
    for row in wanted_bioentries:
        try:
            final_seqfeatures.extend(seqfeature_bioentry_map[row[0]])
        except KeyError:
            pass

    return final_seqfeatures
