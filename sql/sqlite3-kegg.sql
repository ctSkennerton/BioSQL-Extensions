-- the modules come from the kegg database. Kegg does not label their internal
-- nodes in the hierarchy only the leaves so some modules may have an external_id
-- while others (the internal nodes) will not
create table kegg_modules(
    id          integer primary key,
    external_id text,
    name        text
);

-- definition of a closure table for hierarchical data
-- http://www.slideshare.net/billkarwin/models-for-hierarchical-data
create table kegg_module_ontology(
    ancestor   integer not null,
    descendant integer not null,
    length     integer,
    primary key(ancestor, descendant),
    foreign key(ancestor) references kegg_modules(id),
    foreign key(descendant) references kegg_modules(id)
);

create table kegg_pathways(
    id          integer primary key,
    external_id text,
    name        text
);

-- definition of a closure table for hierarchical data
-- http://www.slideshare.net/billkarwin/models-for-hierarchical-data
create table kegg_pathway_ontology(
    ancestor   integer not null,
    descendant integer not null,
    length     integer,
    primary key(ancestor, descendant),
    foreign key(ancestor) references kegg_modules(id),
    foreign key(descendant) references kegg_modules(id)
);

-- a list of kegg orthologs
create table kegg_orthologs(
    id          integer primary key,
    external_id text,
    name        text,
);

create table kegg_module_orthologs(
    ortholog_id integer,
    module_id   integer
);

create table kegg_pathway_orthologs(
    ortholog_id integer,
    module_id   integer
);
