-- BioSQL has a inbuilt ways to add in external database identifiers to a sequence
-- using the dbxref tables. When adding on a KEGG identifier to BioSQL use those
-- inbuilt features to add in the external feature and then use the
-- dbxref_feature_qualifier table to add in extra information about the KEGG entry
-- just as you would a normal seqfeature.

-- These tables are here to help add in the ontology of KEGG pathways and modules
-- and a store of which KEGG orthologs belong to which. I think that BioSQL could
-- do all of this with it's term_relationship and term path tables, but I can't
-- figure out what they are doing. So here I am reinventing the wheel  



-- the modules come from the kegg database. Kegg does not label their internal
-- nodes in the hierarchy only the leaves so some modules may have an external_id
-- while others (the internal nodes) will not
CREATE TABLE kegg_modules(
    id          INTEGER primary key,
    external_id TEXT,
    name        TEXT
);


CREATE TABLE kegg_pathways(
    id          INTEGER primary key,
    external_id TEXT,
    name        TEXT
);

-- definition of a closure table for hierarchical data
-- http://www.slideshare.net/billkarwin/models-for-hierarchical-data
CREATE TABLE kegg_module_ontology(
    ancestor   INTEGER NOT NULL,
    descendant INTEGER NOT NULL,
    length     INTEGER,
    primary key(ancestor, descendant),
    foreign key(ancestor) references kegg_modules(id),
    foreign key(descendant) references kegg_modules(id)
);

 CREATE TABLE kegg_pathway_ontology(
     ancestor   INTEGER NOT NULL,
     descendant INTEGER NOT NULL,
     length     INTEGER,
     primary key(ancestor, descendant),
     foreign key(ancestor) references kegg_pathways(id),
     foreign key(descendant) references kegg_pathways(id)
 );

CREATE TABLE kegg_module_orthologs(
    ortholog_id INTEGER,
    module_id   INTEGER
);

CREATE TABLE kegg_pathway_orthologs(
    ortholog_id INTEGER,
    module_id   INTEGER
);
