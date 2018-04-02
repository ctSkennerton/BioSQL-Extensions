-- get the full taxonomy path for all leaves (species/strains) listed
-- in the taxon and taxon_name tables.
-- creates a table of taxon_id,taxon_string where the taxon_string
-- is 'level1;level2;level3;level4;level5;level6'
CREATE VIEW lineage AS 
SELECT lineage.id,
  string_agg(lineage.name::text, ';'::text) AS lineage
  FROM ( SELECT child.taxon_id AS id, name.name
         FROM taxon child
         JOIN taxon ancestor ON child.left_value >= ancestor.left_value AND child.left_value <= ancestor.right_value
         JOIN taxon_name name ON name.taxon_id = ancestor.taxon_id
         WHERE child.right_value = (child.left_value + 1) AND name.name_class::text = 'scientific name'::text
         ORDER BY ancestor.left_value) AS lineage
  GROUP BY lineage.id;
