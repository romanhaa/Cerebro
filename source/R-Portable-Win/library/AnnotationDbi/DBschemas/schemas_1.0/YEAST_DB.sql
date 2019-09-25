--
-- YEAST_DB schema
-- ===============
--

-- The "sgd" table is the central table.
CREATE TABLE sgd (
  _id INTEGER PRIMARY KEY,
  systematic_name VARCHAR(14) NULL UNIQUE,      -- Yeast ORF ID
  gene_name VARCHAR(14) NULL UNIQUE,            -- Yeast gene name
  sgd_id CHAR(10) NOT NULL UNIQUE               -- SGD ID
);

-- Data linked to the "sgd" table.
CREATE TABLE chromosome_features (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  chromosome VARCHAR(2) NULL,                   -- chromosome name
  start INTEGER NULL,
  feature_description TEXT NULL,                -- Yeast feature description
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE ec (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  ec_number VARCHAR(13) NOT NULL,               -- EC number (no "EC:" prefix)
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE gene2alias (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  alias VARCHAR(13) NOT NULL,                   -- Yeast gene alias
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_bp (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_bp_all (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_cc (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_cc_all (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_mf (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE go_mf_all (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE interpro (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  interpro_id CHAR(9) NOT NULL,                 -- InterPro ID
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE kegg (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  path_id CHAR(5) NOT NULL,                     -- KEGG pathway short ID
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE pfam (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  pfam_id CHAR(7) NOT NULL,                     -- Pfam ID
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE pubmed (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  pubmed_id VARCHAR(10) NOT NULL,               -- PubMed ID
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);
CREATE TABLE smart (
  _id INTEGER NOT NULL,                         -- REFERENCES sgd
  smart_id CHAR(7) NOT NULL,                    -- SMART ID
  FOREIGN KEY (_id) REFERENCES sgd (_id)
);

-- Standalone data tables.
CREATE TABLE chrlengths (
  chromosome VARCHAR(2) PRIMARY KEY,            -- chromosome name
  length INTEGER NOT NULL
);
CREATE TABLE gene2systematic (
  gene_name VARCHAR(14) NULL,                   -- Yeast gene name
  systematic_name VARCHAR(14) NULL              -- Yeast ORF ID
);
CREATE TABLE reject_orf (
  systematic_name VARCHAR(14) PRIMARY KEY       -- Yeast ORF ID
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
CREATE INDEX Fchromosome_features ON chromosome_features (_id);
CREATE INDEX Fec ON ec (_id);
CREATE INDEX Fgene2alias ON gene2alias (_id);
CREATE INDEX Fgo_bp ON go_bp (_id);
CREATE INDEX Fgo_bp_all ON go_bp_all (_id);
CREATE INDEX Fgo_cc ON go_cc (_id);
CREATE INDEX Fgo_cc_all ON go_cc_all (_id);
CREATE INDEX Fgo_mf ON go_mf (_id);
CREATE INDEX Fgo_mf_all ON go_mf_all (_id);
CREATE INDEX Finterpro ON interpro (_id);
CREATE INDEX Fkegg ON kegg (_id);
CREATE INDEX Fpfam ON pfam (_id);
CREATE INDEX Fpubmed ON pubmed (_id);
CREATE INDEX Fsmart ON smart (_id);

