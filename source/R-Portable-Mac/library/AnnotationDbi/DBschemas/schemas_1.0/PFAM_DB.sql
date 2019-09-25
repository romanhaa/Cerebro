CREATE TABLE metadata (
 name TEXT,
 value TEXT
);

CREATE TABLE cazy (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 cazy VARCHAR(6) NOT NULL               --CAZY ID
);

CREATE TABLE homstrad (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 homstrad VARCHAR(20) NOT NULL          --HOMSTRAD ID
);

CREATE TABLE interpro (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 interpro VARCHAR(9) NOT NULL           --INTERPRO ID
);

CREATE TABLE load (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 load VARCHAR(15) NOT NULL              --LOAD ID
);

CREATE TABLE merops (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 merops VARCHAR(3) NOT NULL             --MEROPS ID
);

CREATE TABLE mim (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 mim VARCHAR(6) NOT NULL                --MIM ID
);

CREATE TABLE pdb (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 pdb VARCHAR(6) NOT NULL,            	--PDB ID
 start_point INTEGER,			--start of alignment
 end_point INTEGER			--end of alignment
);

CREATE TABLE pfamb (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 pfamb VARCHAR(8) NOT NULL              --PFAMB ID
);

CREATE TABLE prints (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 prints VARCHAR(7) NOT NULL             --PRINTS ID
);

CREATE TABLE prosite (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 prosite VARCHAR(9) NOT NULL            --PROSITE ID
);

CREATE TABLE prosite_profile (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 prosite_profile VARCHAR(7) NOT NULL    --PROSITE_PROFILE ID
);

CREATE TABLE rm (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 rm VARCHAR(8) NOT NULL                 --RM ID
);

CREATE TABLE scop (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 scop VARCHAR(4) NOT NULL,              --SCOP ID
 placement VARCHAR(2) NOT NULL		--SCOP placement
);

CREATE TABLE smart (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 smart VARCHAR(9) NOT NULL              --SMART ID
);

CREATE TABLE tc (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 tc VARCHAR(6) NOT NULL                 --TC ID
);

CREATE TABLE url (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 url VARCHAR(80) NOT NULL      	        --URL ID
);

CREATE TABLE id (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 id VARCHAR(15) NOT NULL                --PFAM ID ID
);

CREATE TABLE de (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 de VARCHAR(80) NOT NULL                --DE ID
);

CREATE TABLE tp (
 ac VARCHAR(12)  NOT NULL,		--AC ID
 tp VARCHAR(6) NOT NULL                 --TP ID
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
CREATE INDEX cazy_ac ON cazy(ac);
CREATE INDEX homstrad_ac ON homstrad(ac);
CREATE INDEX interpro_ac ON interpro(ac);
CREATE INDEX load_ac ON load(ac);
CREATE INDEX merops_ac ON merops(ac);
CREATE INDEX mim_ac ON mim(ac);
CREATE INDEX pdb_ac ON pdb(ac);
CREATE INDEX pfamb_ac ON pfamb(ac);
CREATE INDEX prints_ac ON prints(ac);
CREATE INDEX prosite_ac ON prosite(ac);
CREATE INDEX prosite_profile_ac ON prosite_profile(ac);
CREATE INDEX rm_ac ON rm(ac);
CREATE INDEX scop_ac ON scop(ac);
CREATE INDEX smart_ac ON smart(ac);
CREATE INDEX tc_ac ON tc(ac);
CREATE INDEX url_ac ON url(ac)
CREATE INDEX id_ac ON id(ac);
CREATE INDEX de_ac ON de(ac);
CREATE INDEX tp_ac ON tp(ac);

