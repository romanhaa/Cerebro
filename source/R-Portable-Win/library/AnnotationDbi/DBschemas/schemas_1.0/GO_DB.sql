--
-- GO_DB schema
-- ============
--

CREATE TABLE go_ontology (
  ontology VARCHAR(9) PRIMARY KEY,              -- GO ontology (short label)
  term_type VARCHAR(18) NOT NULL UNIQUE         -- GO ontology (full label)
);
CREATE TABLE go_obsolete (
  go_id CHAR(10) PRIMARY KEY,                   -- GO ID
  term VARCHAR(255) NOT NULL,                   -- textual label for the GO term
  ontology VARCHAR(9) NOT NULL,                 -- REFERENCES go_ontology
  definition TEXT NULL,                         -- textual definition for the GO term
  FOREIGN KEY (ontology) REFERENCES go_ontology (ontology)
);

-- The "go_term" table is the central table.
CREATE TABLE go_term (
  _id INTEGER PRIMARY KEY, 
  go_id CHAR(10) NOT NULL UNIQUE,               -- GO ID
  term VARCHAR(255) NOT NULL,                   -- textual label for the GO term
  ontology VARCHAR(9) NOT NULL,                 -- REFERENCES go_ontology
  definition TEXT NULL,                         -- textual definition for the GO term
  FOREIGN KEY (ontology) REFERENCES go_ontology (ontology)
);

-- Data linked to the "go_term" table.
CREATE TABLE go_synonym (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  synonym VARCHAR(255) NOT NULL,                -- label or GO ID
  secondary CHAR(10) NULL,                      -- GO ID
  like_go_id SMALLINT,                          -- boolean (1 or 0)
  FOREIGN KEY (_id) REFERENCES go_term (_id)
);

-- child-parent relationship.
CREATE TABLE go_bp_parents (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _parent_id INTEGER NOT NULL,                  -- REFERENCES go_term
  relationship_type VARCHAR(7) NOT NULL,        -- type of GO child-parent relationship
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_parent_id) REFERENCES go_term (_id)
);
CREATE TABLE go_cc_parents (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _parent_id INTEGER NOT NULL,                  -- REFERENCES go_term
  relationship_type VARCHAR(7) NOT NULL,        -- type of GO child-parent relationship
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_parent_id) REFERENCES go_term (_id)
);
CREATE TABLE go_mf_parents (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _parent_id INTEGER NOT NULL,                  -- REFERENCES go_term
  relationship_type VARCHAR(7) NOT NULL,        -- type of GO child-parent relationship
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_parent_id) REFERENCES go_term (_id)
);

-- parent-offspring relationship.
CREATE TABLE go_bp_offspring (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _offspring_id INTEGER NOT NULL,               -- REFERENCES go_term
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_offspring_id) REFERENCES go_term (_id)
);
CREATE TABLE go_cc_offspring (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _offspring_id INTEGER NOT NULL,               -- REFERENCES go_term
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_offspring_id) REFERENCES go_term (_id)
);
CREATE TABLE go_mf_offspring (
  _id INTEGER NOT NULL,                         -- REFERENCES go_term
  _offspring_id INTEGER NOT NULL,               -- REFERENCES go_term
  FOREIGN KEY (_id) REFERENCES go_term (_id),
  FOREIGN KEY (_offspring_id) REFERENCES go_term (_id)
);

-- Metadata tables.
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);
CREATE TABLE map_counts (
  map_name VARCHAR(80) PRIMARY KEY,
  count INTEGER NOT NULL
);
CREATE TABLE map_metadata (
  map_name VARCHAR(80) NOT NULL,
  source_name VARCHAR(80) NOT NULL,
  source_url VARCHAR(255) NOT NULL,
  source_date VARCHAR(20) NOT NULL
);

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX Fgo_term ON go_term (ontology);
CREATE INDEX Fgo_synonym ON go_synonym (_id);
CREATE INDEX Fgo_obsolete ON go_obsolete (ontology);
CREATE INDEX F1go_bp_parents ON go_bp_parents (_id);
CREATE INDEX F2go_bp_parents ON go_bp_parents (_parent_id);
CREATE INDEX F1go_cc_parents ON go_cc_parents (_id);
CREATE INDEX F2go_cc_parents ON go_cc_parents (_parent_id);
CREATE INDEX F1go_mf_parents ON go_mf_parents (_id);
CREATE INDEX F2go_mf_parents ON go_mf_parents (_parent_id);
CREATE INDEX F1go_bp_offspring ON go_bp_offspring (_id);
CREATE INDEX F2go_bp_offspring ON go_bp_offspring (_offspring_id);
CREATE INDEX F1go_cc_offspring ON go_cc_offspring (_id);
CREATE INDEX F2go_cc_offspring ON go_cc_offspring (_offspring_id);
CREATE INDEX F1go_mf_offspring ON go_mf_offspring (_id);
CREATE INDEX F2go_mf_offspring ON go_mf_offspring (_offspring_id);

