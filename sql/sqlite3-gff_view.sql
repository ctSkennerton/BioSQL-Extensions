
-- basic gff information without the attributes text
CREATE VIEW gff_base
AS
    SELECT   e.accession                                  AS fref,
             k2.name                                      AS fsource,
             fl.start_pos                                 AS fstart,
             fl.end_pos                                   AS fend,
             k.name                                       AS type,
             fl.strand                                    AS fstrand,
             f.seqfeature_id                              AS gid
    FROM     seqfeature f,
             location fl,
             bioentry e
    JOIN
             term k
    ON
            f.type_term_id = k.term_id
    JOIN
            term k2
    ON
            f.source_term_id = k2.term_id
    WHERE
             fl.seqfeature_id = f.seqfeature_id           AND
             f.bioentry_id    = e.bioentry_id;

-- this isn't quite a perfect view as the score and phase
-- attributes are left as dummy values. On gff3 files I've loaded
-- these are all stored in the attributes of the record so it
-- would be possible to get them out with a more complicated query
-- but there columns are the least useful bits of information so
-- for now I will leave them be
CREATE VIEW gff
AS
    SELECT   gb.fref                                      AS seqid,
             gb.fsource                                   AS source,
             gb.type                                      AS type,
             gb.fstart                                    AS start,
             gb.fend                                      AS end,
             0                                            AS score,
             gb.fstrand                                   AS strand,
             0                                            AS phase,
             group_concat(t.name || "=" || qv.value, ";") AS attributes
    FROM     gff_base gb
    JOIN     term t,
             seqfeature_qualifier_value qv
    ON       qv.term_id = t.term_id                       AND
             qv.seqfeature_id = gb.gid
    GROUP BY qv.seqfeature_id;
