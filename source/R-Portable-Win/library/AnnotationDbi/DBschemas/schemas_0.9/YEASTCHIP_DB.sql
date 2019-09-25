--
-- YEASTCHIP_DB schema
-- ===================
--

-- The "sgd" table is the central table.
CREATE TABLE sgd (
  id INTEGER PRIMARY KEY,
  systematic_name VARCHAR(14) NULL UNIQUE,      -- Yeast ORF ID
  gene_name VARCHAR(14) NULL UNIQUE,            -- Yeast gene name
  sgd_id CHAR(10) NOT NULL UNIQUE               -- SGD ID
);

-- Data linked to the "sgd" table.
CREATE TABLE probes (
  id INTEGER NULL,                              -- REFERENCES sgd
  probe_id VARCHAR(80) PRIMARY KEY,             -- manufacturer ID
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE chromosome_features (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  chromosome VARCHAR(2) NULL,                   -- chromosome name
  start INTEGER NULL,
  feature_description TEXT NULL,                -- Yeast feature description
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE ec (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  ec_number VARCHAR(13) NOT NULL,               -- EC number (no "EC:" prefix)
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE gene2alias (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  alias VARCHAR(13) NOT NULL,                   -- Yeast gene alias
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_bp (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_bp_all (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_cc (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_cc_all (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_mf (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE go_mf_all (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  go_id CHAR(10) NOT NULL,                      -- GO ID
  evidence CHAR(3) NOT NULL,                    -- GO evidence code
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE kegg (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  kegg_id CHAR(5) NOT NULL,                     -- KEGG pathway short ID
  FOREIGN KEY (id) REFERENCES sgd (id)
);
CREATE TABLE pubmed (
  id INTEGER NOT NULL,                          -- REFERENCES sgd
  pubmed_id VARCHAR(10) NOT NULL,               -- PubMed ID
  FOREIGN KEY (id) REFERENCES sgd (id)
);

-- Standalone data tables.
CREATE TABLE chrlengths (
  chr VARCHAR(2) PRIMARY KEY,                   -- chromosome name
  length INTEGER NOT NULL
);

-- Metadata tables.
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);
CREATE TABLE qcdata (
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
CREATE INDEX Fprobes ON probes (id);
CREATE INDEX Fchromosome_features ON chromosome_features (id);
CREATE INDEX Fec ON ec (id);
CREATE INDEX Fgene2alias ON gene2alias (id);
CREATE INDEX Fgo_bp ON go_bp (id);
CREATE INDEX Fgo_bp_all ON go_bp_all (id);
CREATE INDEX Fgo_cc ON go_cc (id);
CREATE INDEX Fgo_cc_all ON go_cc_all (id);
CREATE INDEX Fgo_mf ON go_mf (id);
CREATE INDEX Fgo_mf_all ON go_mf_all (id);
CREATE INDEX Fkegg ON kegg (id);
CREATE INDEX Fpubmed ON pubmed (id);

