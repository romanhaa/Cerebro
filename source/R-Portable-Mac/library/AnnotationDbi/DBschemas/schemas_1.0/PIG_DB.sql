--
-- PIG_DB schema
-- ===============
--

-- The "genes" table is the central table.
CREATE TABLE genes (
  _id INTEGER PRIMARY KEY,
  gene_id VARCHAR(10) NOT NULL UNIQUE           -- Entrez Gene ID
);

-- Data linked to the "genes" table.
CREATE TABLE accessions (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  accession VARCHAR(20) NOT NULL,               -- GenBank accession number
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE alias (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  alias_symbol VARCHAR(80) NOT NULL,            -- gene symbol or alias
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE chromosomes (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  chromosome VARCHAR(2) NOT NULL,               -- chromosome name
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE chromosome_locations (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  seqname VARCHAR(20) NOT NULL,                 -- sequence name
  start_location INTEGER NOT NULL,
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE ec (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  ec_number VARCHAR(13) NOT NULL,               -- EC number (no "EC:" prefix)
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE gene_info (
  _id INTEGER NOT NULL UNIQUE,                  -- REFERENCES genes
  gene_name VARCHAR(255) NOT NULL,              -- gene name
  symbol VARCHAR(80) NOT NULL,                  -- gene symbol
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE kegg (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  path_id CHAR(5) NOT NULL,                     -- KEGG pathway short ID
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE pubmed (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  pubmed_id VARCHAR(10) NOT NULL,               -- PubMed ID
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE refseq (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  accession VARCHAR(20) NOT NULL,               -- RefSeq accession number
  FOREIGN KEY (_id) REFERENCES genes (_id)
);
CREATE TABLE unigene (
  _id INTEGER NOT NULL,                         -- REFERENCES genes
  unigene_id VARCHAR(10) NOT NULL,              -- UniGene ID
  FOREIGN KEY (_id) REFERENCES genes (_id)
);

-- Standalone data tables.
CREATE TABLE chrlengths (
  chromosome VARCHAR(2) PRIMARY KEY,            -- chromosome name
  length INTEGER NOT NULL
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
CREATE INDEX Faccessions ON accessions (_id);
CREATE INDEX Falias ON alias (_id);
CREATE INDEX Fchromosomes ON chromosomes (_id);
CREATE INDEX Fchromosome_locations ON chromosome_locations (_id);
CREATE INDEX Fec ON ec (_id);
CREATE INDEX Fkegg ON kegg (_id);
CREATE INDEX Fpubmed ON pubmed (_id);
CREATE INDEX Frefseq ON refseq (_id);
CREATE INDEX Funigene ON unigene (_id);

