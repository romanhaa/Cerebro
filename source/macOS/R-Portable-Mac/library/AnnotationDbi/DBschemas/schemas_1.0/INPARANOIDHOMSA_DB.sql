--
-- INPARANOID_DB schema
-- =============
--

CREATE TABLE aedes_aegypti (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE anopheles_gambiae (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE apis_mellifera (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE arabidopsis_thaliana (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE bos_taurus (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE caenorhabditis_briggsae (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE caenorhabditis_elegans (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE caenorhabditis_remanei (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE candida_glabrata (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE canis_familiaris (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE ciona_intestinalis (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE cryptococcus_neoformans (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE danio_rerio (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE debaryomyces_hanseneii (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE dictyostelium_discoideum (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE drosophila_melanogaster (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE drosophila_pseudoobscura (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE entamoeba_histolytica (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE escherichia_coliK12 (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE gallus_gallus (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE gasterosteus_aculeatus (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE kluyveromyces_lactis (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE macaca_mulatta (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE monodelphis_domestica (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE mus_musculus (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE oryza_sativa (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE pan_troglodytes (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE rattus_norvegicus (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE saccharomyces_cerevisiae (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE schizosaccharomyces_pombe (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE takifugu_rubripes (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE tetraodon_nigroviridis (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE xenopus_tropicalis (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);

CREATE TABLE yarrowia_lipolytica (
 inp_id VARCHAR(30) UNIQUE NOT NULL,	--Inparanoid ID
 clust_id INTEGER NOT NULL,		--Inparanoid Cluster ID
 species CHAR(5) NOT NULL,		--Inparanoid Species ID
 score VARCHAR(6) NOT NULL,		--Inparanoid Score
 seed_status CHAR(4)			--Inparanoid Seed Status
);



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
CREATE INDEX aedae_c ON aedes_aegypti(clust_id);
CREATE INDEX aedae_s ON aedes_aegypti(species);
CREATE INDEX anoga_c ON anopheles_gambiae(clust_id);
CREATE INDEX anoga_s ON anopheles_gambiae(species);
CREATE INDEX apime_c ON apis_mellifera(clust_id);
CREATE INDEX apime_s ON apis_mellifera(species);
CREATE INDEX arath_c ON arabidopsis_thaliana(clust_id);
CREATE INDEX arath_s ON arabidopsis_thaliana(species);
CREATE INDEX bosta_c ON bos_taurus(clust_id);
CREATE INDEX bosta_s ON bos_taurus(species);
CREATE INDEX caebr_c ON caenorhabditis_briggsae(clust_id);
CREATE INDEX caebr_s ON caenorhabditis_briggsae(species);
CREATE INDEX caeel_c ON caenorhabditis_elegans(clust_id);
CREATE INDEX caeel_s ON caenorhabditis_elegans(species);
CREATE INDEX caere_c ON caenorhabditis_remanei(clust_id);
CREATE INDEX caere_s ON caenorhabditis_remanei(species);
CREATE INDEX cangl_c ON candida_glabrata(clust_id);
CREATE INDEX cangl_s ON candida_glabrata(species);
CREATE INDEX canfa_c ON canis_familiaris(clust_id);
CREATE INDEX canfa_s ON canis_familiaris(species);
CREATE INDEX cioin_c ON ciona_intestinalis(clust_id);
CREATE INDEX cioin_s ON ciona_intestinalis(species);
CREATE INDEX cryne_c ON cryptococcus_neoformans(clust_id);
CREATE INDEX cryne_s ON cryptococcus_neoformans(species);
CREATE INDEX danre_c ON danio_rerio(clust_id);
CREATE INDEX danre_s ON danio_rerio(species);
CREATE INDEX debha_c ON debaryomyces_hanseneii(clust_id);
CREATE INDEX debha_s ON debaryomyces_hanseneii(species);
CREATE INDEX dicdi_c ON dictyostelium_discoideum(clust_id);
CREATE INDEX dicdi_s ON dictyostelium_discoideum(species);
CREATE INDEX drome_c ON drosophila_melanogaster(clust_id);
CREATE INDEX drome_s ON drosophila_melanogaster(species);
CREATE INDEX drops_c ON drosophila_pseudoobscura(clust_id);
CREATE INDEX drops_s ON drosophila_pseudoobscura(species);
CREATE INDEX enthi_c ON entamoeba_histolytica(clust_id);
CREATE INDEX enthi_s ON entamoeba_histolytica(species);
CREATE INDEX escco_c ON escherichia_coliK12(clust_id);
CREATE INDEX escco_s ON escherichia_coliK12(species);
CREATE INDEX galga_c ON gallus_gallus(clust_id);
CREATE INDEX galga_s ON gallus_gallus(species);
CREATE INDEX gasac_c ON gasterosteus_aculeatus(clust_id);
CREATE INDEX gasac_s ON gasterosteus_aculeatus(species);
CREATE INDEX klula_c ON kluyveromyces_lactis(clust_id);
CREATE INDEX klula_s ON kluyveromyces_lactis(species);
CREATE INDEX macmu_c ON macaca_mulatta(clust_id);
CREATE INDEX macmu_s ON macaca_mulatta(species);
CREATE INDEX mondo_c ON monodelphis_domestica(clust_id);
CREATE INDEX mondo_s ON monodelphis_domestica(species);
CREATE INDEX musmu_c ON mus_musculus(clust_id);
CREATE INDEX musmu_s ON mus_musculus(species);
CREATE INDEX orysa_c ON oryza_sativa(clust_id);
CREATE INDEX orysa_s ON oryza_sativa(species);
CREATE INDEX pantr_c ON pan_troglodytes(clust_id);
CREATE INDEX pantr_s ON pan_troglodytes(species);
CREATE INDEX ratno_c ON rattus_norvegicus(clust_id);
CREATE INDEX ratno_s ON rattus_norvegicus(species);
CREATE INDEX sacce_c ON saccharomyces_cerevisiae(clust_id);
CREATE INDEX sacce_s ON saccharomyces_cerevisiae(species);
CREATE INDEX schpo_c ON schizosaccharomyces_pombe(clust_id);
CREATE INDEX schpo_s ON schizosaccharomyces_pombe(species);
CREATE INDEX fugru_c ON takifugu_rubripes(clust_id);
CREATE INDEX fugru_s ON takifugu_rubripes(species);
CREATE INDEX tetni_c ON tetraodon_nigroviridis(clust_id);
CREATE INDEX tetni_s ON tetraodon_nigroviridis(species);
CREATE INDEX xentr_c ON xenopus_tropicalis(clust_id);
CREATE INDEX xentr_s ON xenopus_tropicalis(species);
CREATE INDEX yarli_c ON yarrowia_lipolytica(clust_id);
CREATE INDEX yarli_s ON yarrowia_lipolytica(species);
