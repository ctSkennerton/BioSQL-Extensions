-- The biosql data model treats a top level record as an individual
-- sequence like a contig as a separate entity. But in the world of
-- draft genomes we usually want to link these together as a genome
-- in many pieces. This is a simple mapping of bioentries that are
-- found in a BioSQL database to their genome


-- Store some basic metadata about our genomes. For now just the name,
-- description and the taxon id should be good enough. Also included
-- is a temporary place to store a taxonomy string if the genome is
-- currently not in genbank and therefore would not have a taxon id.
-- Since the taxon info of BioSQL already stores the name and the tax
-- string in various parts the name, and tmp_tax_string would only
-- needed to be used for novel genomes with the taxon_id field left
-- blank.
CREATE TABLE genome(
    genome_id             INTEGER PRIMARY KEY,
    taxon_id              INTEGER,
    name                  TEXT,
    description           TEXT,
    tmp_tax_string        TEXT
);

-- The BioSQL has a table already called bioentry_relationship which can be
-- used to do something like this but I want another table to keep these
-- higher order mappings separate from the bioentry table
CREATE TABLE genome_sequence(
    genome_id             INTEGER,
    bioentry_id           INTEGER,
    FOREIGN KEY(genome_id) REFERENCES genome(genome_id),
    FOREIGN KEY(bioentry_id) REFERENCES bioentry(bioentry_id)
);
