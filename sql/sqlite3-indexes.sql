-- these are some extra indexes that I've found useful to add onto the default BioSQL tables


-- alot of the time you want to be searching though the values of qualifiers
-- for things like protein IDs, locus_tags that sort of thing. The default
-- index for the seqfeature_qualifier_value uses only the term_id for filtering
-- this index adds the value of the qualifier in as well, which greatly speeds up
-- the process. Of course it will also take up a lot of space for the index
CREATE INDEX seqfeaturequal_type_value ON seqfeature_qualifier_value(term_id, value);
