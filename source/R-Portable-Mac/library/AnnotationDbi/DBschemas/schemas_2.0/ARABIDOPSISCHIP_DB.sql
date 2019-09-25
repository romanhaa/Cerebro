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
CREATE TABLE sqlite_stat1(tbl,idx,stat);
CREATE TABLE probes (probe_id VARCHAR(80) NOT NULL, is_multiple SMALLINT NOT NULL, gene_id VARCHAR(10) NULL);

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX Fprobes ON probes (probe_id);
CREATE INDEX Fgenes ON probes (gene_id);
