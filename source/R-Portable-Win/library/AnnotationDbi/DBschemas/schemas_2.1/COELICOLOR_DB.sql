CREATE TABLE gene_info (
 _id INTEGER REFERENCES genes(_id),
 gene_name TEXT,
 symbol TEXT
);
CREATE TABLE refseq (
 _id INTEGER REFERENCES genes(_id),
 accession TEXT
);
CREATE TABLE ec (
 _id INTEGER REFERENCES genes(_id),
 ec_number TEXT
);
CREATE TABLE kegg (
 _id INTEGER REFERENCES genes(_id),
 path_id TEXT
);
CREATE TABLE map_metadata (
 map_name TEXT,
 source_name TEXT,
 source_url TEXT,
 source_date TEXT
);
CREATE TABLE pubmed (
 _id INTEGER REFERENCES genes(_id),
 pubmed_id TEXT
);
CREATE TABLE locus (
 _id INTEGER REFERENCES genes(_id),
 locus_tag TEXT                            -- SCO IDs
);
CREATE TABLE genes (
 _id INTEGER PRIMARY KEY,
 gene_id TEXT                            -- Entrez_Gene IDs
);
CREATE TABLE accessions (
 _id INTEGER NOT NULL,                     
 accession TEXT,                           -- RefSeq accessions
 FOREIGN KEY (_id) REFERENCES genes(_id)
);
CREATE TABLE protein (
 _id INTEGER REFERENCES genes(_id),
 protein_gi TEXT
);
CREATE TABLE chrlengths (
      chromosome VARCHAR(1) PRIMARY KEY,            -- chromosome name [0: chromosome, 1: plasmid SCP1, 2: plasmid SCP2]
      length INTEGER NOT NULL
    );
CREATE TABLE chromosomes (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      chromosome VARCHAR(1),                        -- chromosome name
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE alias (
      _id INTEGER NOT NULL,                         -- REFERENCES  genes 
      alias_symbol VARCHAR(80) NOT NULL,            -- gene symbol or alias (synonyms)
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
);
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);
CREATE TABLE go_bp (
 _id INTEGER REFERENCES genes(_id),
 go_id TEXT,
 evidence TEXT
);
CREATE TABLE go_cc (
 _id INTEGER REFERENCES genes(_id),
 go_id TEXT,
 evidence TEXT
);
CREATE TABLE go_mf (
 _id INTEGER REFERENCES genes(_id),
 go_id TEXT,
 evidence TEXT
);
CREATE TABLE go_bp_all (_id INTEGER REFERENCES genes(_id),go_id TEXT,evidence TEXT);
CREATE TABLE go_cc_all (_id INTEGER REFERENCES genes(_id),go_id TEXT,evidence TEXT);
CREATE TABLE go_mf_all (_id INTEGER REFERENCES genes(_id),go_id TEXT,evidence TEXT);
CREATE TABLE map_counts (
      map_name VARCHAR(80) PRIMARY KEY,
      count INTEGER NOT NULL
    );
CREATE TABLE chromosome_locations (
      _id INTEGER NOT NULL,                      -- REFERENCES  genes
      seqname VARCHAR(20) NOT NULL,           -- chromosome name
      start_location INTEGER NOT NULL,
      end_location INTEGER NOT NULL,
      FOREIGN KEY (_id) REFERENCES  genes  (_id)
    );
CREATE TABLE symbol (
 _id INTEGER REFERENCES genes(_id),
 symbol TEXT
);

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX en1 ON ec(_id);
CREATE INDEX gi1 ON gene_info(_id);
CREATE INDEX k1 ON kegg(_id);
CREATE INDEX pu1 ON pubmed(_id);
CREATE INDEX rs1 on refseq(_id);
CREATE INDEX lt1 ON locus(_id);
CREATE INDEX lt2 ON locus(locus_tag);
CREATE INDEX ge2 ON genes(gene_id);
CREATE INDEX ge1 on genes(_id);
CREATE INDEX ac1 ON accessions(_id);
CREATE INDEX ac2 ON accessions(accession);
CREATE INDEX pr1 ON protein(_id);
CREATE INDEX al1 ON alias(_id);
CREATE INDEX ch ON chromosomes (_id);
CREATE INDEX go1 ON go_bp(_id);
CREATE INDEX go2 ON go_mf(_id);
CREATE INDEX go3 ON go_cc(_id);
CREATE INDEX go4 ON go_bp(go_id);
CREATE INDEX go5 ON go_mf(go_id);
CREATE INDEX go6 ON go_cc(go_id);
CREATE INDEX chl ON chromosome_locations (_id);
CREATE INDEX gisymbol ON symbol(_id);
