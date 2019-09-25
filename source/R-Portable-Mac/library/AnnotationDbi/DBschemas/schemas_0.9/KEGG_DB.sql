--
-- KEGG_DB schema
-- ==============
--

CREATE TABLE ec2go (
  ec_no VARCHAR(16) NOT NULL,                   -- EC number (with "EC:" prefix)
  go_id CHAR(10) NOT NULL                       -- GO ID
);
CREATE TABLE pathway2gene (
  pathway_id CHAR(8) NOT NULL,                  -- KEGG pathway long ID
  gene_id VARCHAR(20) NOT NULL                  -- Entrez Gene or ORF ID
);
CREATE TABLE pathway2name (
  path_id CHAR(5) NOT NULL UNIQUE,              -- KEGG pathway short ID 
  path_name VARCHAR(80) NOT NULL UNIQUE         -- KEGG pathway name
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

-- Indexes.
CREATE INDEX Ipathway2gene ON pathway2gene (gene_id);

