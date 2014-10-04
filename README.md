This is a repository of sql tables and views that could be added onto an
existing BioSQL database. They are written only for the sqlite database engine
and are primarily related to adding functionality for performing proteomics
experiments. I have not tried to write much of this in a portable way; the
tables were created for my particular purposes. So the restrictions are that you
are using Sipros for your proteome searches, cause that is what we use. There
is some more general stuff like a couple of tables for higher order mapping of
bioentries into genomes, which is useful for when you have a draft genome
broken into many contigs (which each map to a single bioentry). There is also
a semi-correct gff3 view that could be added to any sqlite3 BioSQL database.
