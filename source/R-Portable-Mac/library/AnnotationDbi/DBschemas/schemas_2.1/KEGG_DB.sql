CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);
CREATE TABLE pathway2name (
  path_id CHAR(5) NOT NULL UNIQUE,              -- KEGG pathway short ID
  path_name VARCHAR(100) NOT NULL UNIQUE         -- KEGG pathway name
);
CREATE TABLE pathway2gene (
  pathway_id CHAR(8) NOT NULL,                  -- KEGG pathway long ID
  gene_or_orf_id VARCHAR(20) NOT NULL                  -- Entrez Gene or ORF ID
);
CREATE TABLE ath_NCBI_pathway2gene (
  pathway_id CHAR(8) NOT NULL,                  -- KEGG pathway long ID
  gene_or_orf_id VARCHAR(20) NOT NULL                  -- Entrez Gene or ORF ID
);
CREATE TABLE ec2go(ec_no VARCHAR(16),go_id CHAR(10));
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
CREATE INDEX Ipathway2gene ON pathway2gene (gene_or_orf_id);
CREATE INDEX ath_NCBI_Ipathway2gene ON ath_NCBI_pathway2gene (gene_or_orf_id);
