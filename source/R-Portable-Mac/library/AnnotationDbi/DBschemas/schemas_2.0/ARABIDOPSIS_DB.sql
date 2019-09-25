CREATE TABLE metadata (name VARCHAR(80) PRIMARY KEY, value VARCHAR(255) );
CREATE TABLE map_metadata (
      map_name VARCHAR(80) NOT NULL,
      source_name VARCHAR(80) NOT NULL,
      source_url VARCHAR(255) NOT NULL,
      source_date VARCHAR(20) NOT NULL
    );
CREATE TABLE map_counts (
      map_name VARCHAR(80) PRIMARY KEY,
      count INTEGER NOT NULL
    );
CREATE TABLE genes (
      _id INTEGER PRIMARY KEY,
      gene_id CHAR(9) NOT NULL UNIQUE               -- AGI locus ID
    );
CREATE TABLE sqlite_stat1(tbl,idx,stat);
CREATE TABLE entrez_genes (
      _id INTEGER PRIMARY KEY,
      gene_id CHAR(9) NOT NULL UNIQUE               -- AGI locus ID
    );
CREATE TABLE gene_info (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      gene_name VARCHAR(255) NULL,                  -- gene name
      symbol VARCHAR(80) NULL,                      -- gene symbol
      chromosome CHAR(1) NULL,                      -- Arabidopsis chromosome
      FOREIGN KEY (_id) REFERENCES  genes  (_id));
CREATE TABLE pubmed (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      pubmed_id VARCHAR(10) NOT NULL,               -- PubMed ID
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE chrlengths (
      chromosome VARCHAR(2) PRIMARY KEY,                   -- chromosome name
      length INTEGER NOT NULL
    );
CREATE TABLE go_bp (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE go_mf (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE go_cc (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE go_bp_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE go_mf_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE go_cc_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE kegg (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      path_id CHAR(5) NOT NULL,                     -- KEGG pathway short ID
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE ec (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      ec_number VARCHAR(13) NOT NULL,               -- EC number 
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE chromosome_locations (
      _id INTEGER NOT NULL,                      -- REFERENCES  genes 
      seqname VARCHAR(20) NOT NULL,              -- sequence name
      start_location INTEGER NOT NULL,
      end_location INTEGER NOT NULL,
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE refseq (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      accession VARCHAR(20) NOT NULL,               -- RefSeq accession number
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE aracyc (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      pathway_name VARCHAR(255) NOT NULL,           -- AraCyc pathway name
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE enzyme (
      _id INTEGER NOT NULL,                          -- REFERENCES  genes 
      ec_name VARCHAR(255) NOT NULL,                -- EC name
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX Fentrez_genes ON entrez_genes(gene_id);
CREATE INDEX Fgene_info ON gene_info (_id);
CREATE INDEX Fpubmed ON pubmed (_id);
CREATE INDEX Fgo_bp ON go_bp (_id);
CREATE INDEX Fgo_bp_go_id ON go_bp (go_id);
CREATE INDEX Fgo_mf ON go_mf (_id);
CREATE INDEX Fgo_mf_go_id ON go_mf (go_id);
CREATE INDEX Fgo_cc ON go_cc (_id);
CREATE INDEX Fgo_cc_go_id ON go_cc (go_id);
CREATE INDEX Fgo_bp_all ON go_bp_all (_id);
CREATE INDEX Fgo_bp_all_go_id ON go_bp_all (go_id);
CREATE INDEX Fgo_mf_all ON go_mf_all (_id);
CREATE INDEX Fgo_mf_all_go_id ON go_mf_all (go_id);
CREATE INDEX Fgo_cc_all ON go_cc_all (_id);
CREATE INDEX Fgo_cc_all_go_id ON go_cc_all (go_id);
CREATE INDEX Fkegg ON kegg (_id);
CREATE INDEX Fec ON ec (_id);
CREATE INDEX Fchromosome_locations ON chromosome_locations (_id);
CREATE INDEX Frefseq ON refseq (_id);
CREATE INDEX Faracyc ON aracyc (_id);
CREATE INDEX Fenzyme ON enzyme (_id);
