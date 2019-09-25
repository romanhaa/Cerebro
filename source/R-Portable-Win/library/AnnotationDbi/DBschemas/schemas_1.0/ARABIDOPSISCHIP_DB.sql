--
-- ARABIDOPSISCHIP_DB schema
-- =========================
--

-- The "genes" table is the central table.
CREATE TABLE genes (
  _id INTEGER PRIMARY KEY,
  gene_id CHAR(9) NOT NULL UNIQUE               -- AGI locus ID
);

-- Data linked to the "genes" table.
CREATE TABLE probes (
  probe_id VARCHAR(80) NOT NULL,                -- manufacturer ID
  is_multiple SMALLINT NOT NULL,                -- a silly and useless field
  _id INTEGER NULL,                             -- REFERENCES genes
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE aracyc (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  pathway_name VARCHAR(255) NOT NULL,            -- AraCyc pathway name
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE chromosome_locations (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  seqname CHAR(1) NOT NULL,                      -- Arabidopsis chromosome
  start_location INTEGER NOT NULL,
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE ec (                                --    Table
  _id INTEGER NOT NULL,                          -- 
  ec_number VARCHAR(13) NOT NULL,                --    NOT
  FOREIGN KEY (_id) REFERENCES genes (_id)       -- 
);                                               --    used!
CREATE TABLE enzyme (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  ec_name VARCHAR(255) NOT NULL,                 -- EC name
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
-- Note that the "gene_info" table differs from other schemas:
--   o no UNIQUE constraint on col "id"
--   o no NOT NULL constraints on cols "gene_name" and "symbol"
--   o one additional col "chromosome"
CREATE TABLE gene_info (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  gene_name VARCHAR(255) NULL,                   -- gene name
  symbol VARCHAR(80) NULL,                       -- gene symbol
  chromosome CHAR(1) NULL,                       -- Arabidopsis chromosome
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_bp (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_bp_all (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_cc (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_cc_all (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_mf (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE go_mf_all (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  go_id CHAR(10) NOT NULL,                       -- GO ID
  evidence CHAR(3) NOT NULL,                     -- GO evidence code
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE kegg (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  path_id CHAR(5) NOT NULL,                      -- KEGG pathway short ID
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE pubmed (
  _id INTEGER NOT NULL,                          -- REFERENCES genes
  pubmed_id VARCHAR(10) NOT NULL,                -- PubMed ID
  FOREIGN KEY (_id) REFERENCES genes (_id)
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
CREATE INDEX Fprobes ON probes (_id);
CREATE INDEX Faracyc ON aracyc (_id);
CREATE INDEX Fchromosome_locations ON chromosome_locations (_id);
CREATE INDEX Fec ON ec (_id);
CREATE INDEX Fenzyme ON enzyme (_id);
CREATE INDEX Fgene_info ON gene_info (_id);
CREATE INDEX Fgo_bp ON go_bp (_id);
CREATE INDEX Fgo_bp_all ON go_bp_all (_id);
CREATE INDEX Fgo_cc ON go_cc (_id);
CREATE INDEX Fgo_cc_all ON go_cc_all (_id);
CREATE INDEX Fgo_mf ON go_mf (_id);
CREATE INDEX Fgo_mf_all ON go_mf_all (_id);
CREATE INDEX Fkegg ON kegg (_id);
CREATE INDEX Fpubmed ON pubmed (_id);

-- Other indexes.
CREATE INDEX Lprobes ON probes (probe_id);

