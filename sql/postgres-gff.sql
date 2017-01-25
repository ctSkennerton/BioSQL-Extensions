SELECT seqfeature.bioentry_id,
    k2.name AS source,
    location.start_pos,
    location.end_pos,
    k1.name,
    NULL::unknown AS score,
    location.strand,
    NULL::unknown AS phase,
    seqfeature.seqfeature_id,
    string_agg((term.name::text || '='::text) || seqfeature_qualifier_value.value, ';'::text ORDER BY term.name) AS attributes
   FROM seqfeature
     JOIN term k1 ON k1.term_id = seqfeature.type_term_id
     JOIN term k2 ON k2.term_id = seqfeature.source_term_id
     JOIN location ON location.seqfeature_id = seqfeature.seqfeature_id
     JOIN seqfeature_qualifier_value ON seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
     JOIN term ON seqfeature_qualifier_value.term_id = term.term_id AND term.name::text <> 'source'::text
  GROUP BY seqfeature.bioentry_id, k2.name, location.start_pos, location.end_pos, k1.name, location.strand, seqfeature.seqfeature_id, seqfeature_qualifier_value.seqfeature_id; 
