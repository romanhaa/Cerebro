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
CREATE TABLE sgd (
      _id INTEGER PRIMARY KEY,
      systematic_name VARCHAR(14) NULL UNIQUE,      -- Yeast gene systematic name
      gene_name VARCHAR(14) NULL UNIQUE,            -- Yeast gene name
      sgd_id CHAR(10) NOT NULL UNIQUE               -- SGD ID
    );
CREATE TABLE sqlite_stat1(tbl,idx,stat);
CREATE TABLE chromosome_features (
        _id INTEGER NOT NULL,                         -- REFERENCES sgd
        chromosome VARCHAR(2) NULL,                   -- chromosome name
        start INTEGER NULL,
        stop INTEGER NULL,
        feature_description TEXT NULL,                -- Yeast feature description
        FOREIGN KEY (_id) REFERENCES sgd (_id)
      );
CREATE TABLE gene2alias (
      _id INTEGER NOT NULL,                         -- REFERENCES sgd
      alias VARCHAR(13) NOT NULL,                   -- Yeast gene alias
      FOREIGN KEY (_id) REFERENCES sgd (_id)
    );
CREATE TABLE pubmed (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      pubmed_id VARCHAR(10) NOT NULL,               -- PubMed ID
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE chrlengths (
      chromosome VARCHAR(2) PRIMARY KEY,                   -- chromosome name
      length INTEGER NOT NULL
    );
CREATE TABLE go_bp (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE go_mf (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE go_cc (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE go_bp_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE go_mf_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE go_cc_all (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      go_id CHAR(10) NOT NULL,                      -- GO ID
      evidence CHAR(3) NOT NULL,                    -- GO evidence code
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE kegg (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      path_id CHAR(5) NOT NULL,                     -- KEGG pathway short ID
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ec (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      ec_number VARCHAR(13) NOT NULL,               -- EC number 
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE pfam (
      _id INTEGER NOT NULL,                         -- REFERENCES sgd
      pfam_id CHAR(7) NOT NULL,                     -- Pfam ID
      FOREIGN KEY (_id) REFERENCES sgd (_id));
CREATE TABLE smart (
      _id INTEGER NOT NULL,                         -- REFERENCES sgd 
      smart_id CHAR(7) NOT NULL,                    -- SMART ID
      FOREIGN KEY (_id) REFERENCES sgd (_id));
CREATE TABLE interpro (
      _id INTEGER NOT NULL,                         -- REFERENCES sgd
      interpro_id CHAR(9) NOT NULL,                 -- InterPro ID
      FOREIGN KEY (_id) REFERENCES sgd (_id));
CREATE TABLE reject_orf (
      systematic_name VARCHAR(14) PRIMARY KEY);
CREATE TABLE gene2systematic (
      gene_name VARCHAR(14) NULL,                     -- Yeast gene name
      systematic_name VARCHAR(14) NULL);
CREATE TABLE genes (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      gene_id VARCHAR(20) NOT NULL,                 -- gene id
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ensembl (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      ensembl_id VARCHAR(20) NOT NULL,              -- ensembl id
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ensembl2ncbi (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      ensembl_id VARCHAR(20) NOT NULL,              -- ensembl id
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ncbi2ensembl (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      ensembl_id VARCHAR(20) NOT NULL,              -- ensembl id
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ensembl_prot (
      _id INTEGER NOT NULL,                          -- REFERENCES  sgd 
      prot_id VARCHAR(20) NOT NULL,                  -- Ensembl Protein ID
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE ensembl_trans (
      _id INTEGER NOT NULL,                          -- REFERENCES  sgd 
      trans_id VARCHAR(20) NOT NULL,                  -- Ensembl Transcript ID
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE uniprot (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      uniprot_id VARCHAR(20) NOT NULL,              -- uniprot id
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );
CREATE TABLE refseq (
      _id INTEGER NOT NULL,                         -- REFERENCES  sgd 
      accession VARCHAR(20) NOT NULL,               -- RefSeq accession number
      FOREIGN KEY (_id) REFERENCES  sgd  (_id)
    );

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX Fchromosome_features ON chromosome_features (_id);
CREATE INDEX Fgene2alias ON gene2alias(_id);
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
CREATE INDEX Fpfam ON pfam (_id);
CREATE INDEX Fsmart ON smart (_id);
CREATE INDEX Finterpro ON interpro (_id);
CREATE INDEX Fgene ON genes(_id);
CREATE INDEX Fensembl ON ensembl (_id);
CREATE INDEX Fensembl2ncbi ON ensembl2ncbi (_id);
CREATE INDEX Fncbi2ensembl ON ncbi2ensembl (_id);
CREATE INDEX Fensemblp ON ensembl_prot (_id);
CREATE INDEX Fensemblt ON ensembl_trans (_id);
CREATE INDEX Funiprot ON uniprot (_id);
CREATE INDEX Frefseq ON refseq (_id);
