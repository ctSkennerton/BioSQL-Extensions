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


-- a list of kegg orthologs
CREATE TABLE kegg_orthologs(
    id          INTEGER primary key,
    external_id TEXT,
    gene_id     TEXT,
    product     TEXT,
    ec          TEXT
);

CREATE TABLE kegg_module_orthologs(
    ortholog_id INTEGER,
    module_id   INTEGER
);

CREATE TABLE kegg_pathway_orthologs(
    ortholog_id INTEGER,
    module_id   INTEGER
);
