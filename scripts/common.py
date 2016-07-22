from Bio.Seq import reverse_complement, translate as bio_translate
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

def standard_options():
    import argparse
    from getpass import getuser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', help='name of premade biosql database')
    parser.add_argument('-r', '--driver', help='Python database driver to use (must be installed separately)', choices=["MySQLdb", "psycopg2", "sqlite3"], default='psycopg2')
    parser.add_argument('-p', '--port', help='post to connect to on the host',
            default=5432)
    parser.add_argument('-u', '--user', help='database user name',
            default=getuser())
    parser.add_argument('-P','--password', help='database password for user')
    parser.add_argument('-H', '--host', help='host to connect to',
            default='localhost')

    return parser

def get_seqfeature_id_from_qv(db, qualifier, value, biodatabase_id=None):
    sql = r'select seqfeature_id from seqfeature_qualifier_value join term using(term_id) join seqfeature using(seqfeature_id) join bioentry using(bioentry_id) where term.name = %s and value = %s'
    if biodatabase_id is not None:
        sql += ' and biodatabase_id = %s'

    if biodatabase_id is not None:
        rows = db.adaptor.execute_and_fetchall(sql, (qualifier, value, db.dbid))
    else:
        rows = db.adaptor.execute_and_fetchall(sql, (qualifier, value))

    if len(rows) > 1:
        raise ValueError("There is more than one seqfeature associated with qualifier={} and value={}".format(qualifier, value))
    elif len(rows) == 0:
        raise ValueError("There are no seqfeature associated with qualifier={} and value={}".format(qualifier, value))

    return rows[0][0]


def get_seqfeature_ids_from_qv(db, qualifier, value, biodatabase_id=None):
    sql = r'select seqfeature_id from seqfeature_qualifier_value join term using(term_id) join seqfeature using(seqfeature_id) join bioentry using(bioentry_id) where term.name = %s and value = %s'
    if biodatabase_id is not None:
        sql += ' and biodatabase_id = %s'

    if biodatabase_id is not None:
        col0 = db.adaptor.execute_and_fetch_col0(sql, (qualifier, value, db.dbid))
    else:
        col0 = db.adaptor.execute_and_fetch_col0(sql, (qualifier, value))

    if len(col0) == 0:
        raise ValueError("There are no seqfeature associated with qualifier={} and value={}".format(qualifier, value))

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
    bioentry_ids = []
    for c in chunks(ids, 900):
        sql = "SELECT d.name, s.bioentry_id, s.seqfeature_id FROM seqfeature s  \
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

def extract_feature_sql(server, seqfeature_ids, type=['CDS', 'rRNA', 'tRNA'], qualifier=['ID','locus_tag'], translate=False):
    """raw sql extraction of fasta seqfeatures
    """
    for chunk in chunks(seqfeature_ids, 900):
        sql = "SELECT f.seqfeature_id AS gid, \
                    fl.strand,\
                    substring(s.seq, fl.start_pos, (fl.end_pos - fl.start_pos)+1) AS subseq\
               FROM seqfeature f,\
                    location fl,\
                    biosequence s\
               WHERE f.seqfeature_id IN ({}) AND\
                     fl.seqfeature_id = f.seqfeature_id AND\
                     f.bioentry_id = s.bioentry_id".format(generate_placeholders(len(chunk)))

        features = server.adaptor.execute_and_fetchall(sql, tuple(chunk) )
        results = {}
        for sfid, strand, seq in features:
            results[sfid] = (strand, seq)

        qual_select_sql = 'SELECT seqfeature_id, name, value FROM seqfeature_qualifier_value qv, term t WHERE seqfeature_id IN ({}) AND t.term_id = qv.term_id'.format(generate_placeholders(len(chunk)))
        qv = {}
        for seqfeature_id, name, value in server.adaptor.execute_and_fetchall(qual_select_sql, tuple(chunk)):
            try:
                qv[seqfeature_id][name] = value
            except KeyError:
                qv[seqfeature_id] = {}
                qv[seqfeature_id][name] = value

        for seqfeature_id, (strand, seq) in results.items():
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

            if translate:
                seq = bio_translate(seq)
            try:
                name += ' ' + qv[seqfeature_id]['product']
            except KeyError:
                pass

            print(">{}\n{}".format(name, seq))
