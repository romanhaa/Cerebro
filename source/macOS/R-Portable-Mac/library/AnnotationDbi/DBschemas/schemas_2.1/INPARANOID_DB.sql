CREATE TABLE metadata (
      name VARCHAR(80) PRIMARY KEY,
      value VARCHAR(255));
CREATE TABLE map_metadata (
      map_name VARCHAR(80) NOT NULL,
      source_name VARCHAR(80) NOT NULL,
      source_url VARCHAR(255) NOT NULL,
      source_date VARCHAR(20) NOT NULL);
CREATE TABLE map_counts (
      map_name VARCHAR(80) PRIMARY KEY,
      count INTEGER NOT NULL);
CREATE TABLE Acyrthosiphon_pisum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Aedes_aegypti (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Anopheles_gambiae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Apis_mellifera (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Arabidopsis_thaliana (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Aspergillus_fumigatus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Batrachochytrium_dendrobatidis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Bombyx_mori (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Bos_taurus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Branchiostoma_floridae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Brugia_malayi (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Caenorhabditis_brenneri (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Caenorhabditis_briggsae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Caenorhabditis_elegans (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Caenorhabditis_japonica (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Caenorhabditis_remanei (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Candida_albicans (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Candida_glabrata (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Canis_familiaris (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Capitella_spI (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Cavia_porcellus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Chlamydomonas_reinhardtii (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ciona_intestinalis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ciona_savignyi (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Coccidioides_immitis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Coprinopsis_cinereus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Cryptococcus_neoformans (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Cryptosporidium_hominis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Cryptosporidium_parvum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Culex_pipiens (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Cyanidioschyzon_merolae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Danio_rerio (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Daphnia_pulex (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Debaryomyces_hansenii (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Dictyostelium_discoideum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_ananassae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_grimshawi (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_melanogaster (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_mojavensis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_pseudoobscura (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_virilis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Drosophila_willistoni (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Entamoeba_histolytica (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Equus_caballus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Escherichia_coliK12 (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Fusarium_graminearum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Gallus_gallus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Gasterosteus_aculeatus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Giardia_lamblia (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Helobdella_robusta (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ixodes_scapularis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Kluyveromyces_lactis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Leishmania_major (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Lottia_gigantea (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Macaca_mulatta (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Magnaporthe_grisea (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Monodelphis_domestica (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Monosiga_brevicollis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Mus_musculus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Nasonia_vitripennis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Nematostella_vectensis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Neurospora_crassa (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ornithorhynchus_anatinus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Oryza_sativa (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Oryzias_latipes (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ostreococcus_tauri (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Pan_troglodytes (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Pediculus_humanus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Physcomitrella_patens (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Phytophthora_ramorum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Phytophthora_sojae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Plasmodium_falciparum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Plasmodium_vivax (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Pongo_pygmaeus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Populus_trichocarpa (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Pristionchus_pacificus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Puccinia_graminis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Rattus_norvegicus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Rhizopus_oryzae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Saccharomyces_cerevisiae (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Schistosoma_mansoni (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Schizosaccharomyces_pombe (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Sclerotinia_sclerotiorum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Sorghum_bicolor (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Stagonospora_nodorum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Strongylocentrotus_purpuratus (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Takifugu_rubripes (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Tetrahymena_thermophila (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Tetraodon_nigroviridis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Thalassiosira_pseudonana (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Theileria_annulata (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Theileria_parva (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Tribolium_castaneum (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Trichomonas_vaginalis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Trichoplax_adhaerens (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Trypanosoma_cruzi (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Ustilago_maydis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Xenopus_tropicalis (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));
CREATE TABLE Yarrowia_lipolytica (
        inp_id VARCHAR(30) UNIQUE NOT NULL,
        clust_id INTEGER NOT NULL,
        species CHAR(5) NOT NULL,
        score VARCHAR(6) NOT NULL,
        seed_status CHAR(4));

-- Explicit index creation on the referencing column of all the foreign keys.
-- Note that this is only needed for SQLite: PostgreSQL and MySQL create those
-- indexes automatically.
CREATE INDEX ACYPI_c ON Acyrthosiphon_pisum(clust_id);
CREATE INDEX ACYPI_s ON Acyrthosiphon_pisum(species);
CREATE INDEX AEDAE_c ON Aedes_aegypti(clust_id);
CREATE INDEX AEDAE_s ON Aedes_aegypti(species);
CREATE INDEX ANOGA_c ON Anopheles_gambiae(clust_id);
CREATE INDEX ANOGA_s ON Anopheles_gambiae(species);
CREATE INDEX APIME_c ON Apis_mellifera(clust_id);
CREATE INDEX APIME_s ON Apis_mellifera(species);
CREATE INDEX ARATH_c ON Arabidopsis_thaliana(clust_id);
CREATE INDEX ARATH_s ON Arabidopsis_thaliana(species);
CREATE INDEX ASPFU_c ON Aspergillus_fumigatus(clust_id);
CREATE INDEX ASPFU_s ON Aspergillus_fumigatus(species);
CREATE INDEX BATDE_c ON Batrachochytrium_dendrobatidis(clust_id);
CREATE INDEX BATDE_s ON Batrachochytrium_dendrobatidis(species);
CREATE INDEX BOMMO_c ON Bombyx_mori(clust_id);
CREATE INDEX BOMMO_s ON Bombyx_mori(species);
CREATE INDEX BOSTA_c ON Bos_taurus(clust_id);
CREATE INDEX BOSTA_s ON Bos_taurus(species);
CREATE INDEX BRAFL_c ON Branchiostoma_floridae(clust_id);
CREATE INDEX BRAFL_s ON Branchiostoma_floridae(species);
CREATE INDEX BRUMA_c ON Brugia_malayi(clust_id);
CREATE INDEX BRUMA_s ON Brugia_malayi(species);
CREATE INDEX CAEBRE_c ON Caenorhabditis_brenneri(clust_id);
CREATE INDEX CAEBRE_s ON Caenorhabditis_brenneri(species);
CREATE INDEX CAEBR_c ON Caenorhabditis_briggsae(clust_id);
CREATE INDEX CAEBR_s ON Caenorhabditis_briggsae(species);
CREATE INDEX CAEEL_c ON Caenorhabditis_elegans(clust_id);
CREATE INDEX CAEEL_s ON Caenorhabditis_elegans(species);
CREATE INDEX CAEJA_c ON Caenorhabditis_japonica(clust_id);
CREATE INDEX CAEJA_s ON Caenorhabditis_japonica(species);
CREATE INDEX CAERE_c ON Caenorhabditis_remanei(clust_id);
CREATE INDEX CAERE_s ON Caenorhabditis_remanei(species);
CREATE INDEX CANAL_c ON Candida_albicans(clust_id);
CREATE INDEX CANAL_s ON Candida_albicans(species);
CREATE INDEX CANGL_c ON Candida_glabrata(clust_id);
CREATE INDEX CANGL_s ON Candida_glabrata(species);
CREATE INDEX CANFA_c ON Canis_familiaris(clust_id);
CREATE INDEX CANFA_s ON Canis_familiaris(species);
CREATE INDEX CAPSP_c ON Capitella_spI(clust_id);
CREATE INDEX CAPSP_s ON Capitella_spI(species);
CREATE INDEX CAVPO_c ON Cavia_porcellus(clust_id);
CREATE INDEX CAVPO_s ON Cavia_porcellus(species);
CREATE INDEX CHLRE_c ON Chlamydomonas_reinhardtii(clust_id);
CREATE INDEX CHLRE_s ON Chlamydomonas_reinhardtii(species);
CREATE INDEX CIOIN_c ON Ciona_intestinalis(clust_id);
CREATE INDEX CIOIN_s ON Ciona_intestinalis(species);
CREATE INDEX CIOSA_c ON Ciona_savignyi(clust_id);
CREATE INDEX CIOSA_s ON Ciona_savignyi(species);
CREATE INDEX COCIM_c ON Coccidioides_immitis(clust_id);
CREATE INDEX COCIM_s ON Coccidioides_immitis(species);
CREATE INDEX COPCI_c ON Coprinopsis_cinereus(clust_id);
CREATE INDEX COPCI_s ON Coprinopsis_cinereus(species);
CREATE INDEX CRYNE_c ON Cryptococcus_neoformans(clust_id);
CREATE INDEX CRYNE_s ON Cryptococcus_neoformans(species);
CREATE INDEX CRYHO_c ON Cryptosporidium_hominis(clust_id);
CREATE INDEX CRYHO_s ON Cryptosporidium_hominis(species);
CREATE INDEX CRYPA_c ON Cryptosporidium_parvum(clust_id);
CREATE INDEX CRYPA_s ON Cryptosporidium_parvum(species);
CREATE INDEX CULPI_c ON Culex_pipiens(clust_id);
CREATE INDEX CULPI_s ON Culex_pipiens(species);
CREATE INDEX CYAME_c ON Cyanidioschyzon_merolae(clust_id);
CREATE INDEX CYAME_s ON Cyanidioschyzon_merolae(species);
CREATE INDEX DANRE_c ON Danio_rerio(clust_id);
CREATE INDEX DANRE_s ON Danio_rerio(species);
CREATE INDEX DAPPU_c ON Daphnia_pulex(clust_id);
CREATE INDEX DAPPU_s ON Daphnia_pulex(species);
CREATE INDEX DEBHA_c ON Debaryomyces_hansenii(clust_id);
CREATE INDEX DEBHA_s ON Debaryomyces_hansenii(species);
CREATE INDEX DICDI_c ON Dictyostelium_discoideum(clust_id);
CREATE INDEX DICDI_s ON Dictyostelium_discoideum(species);
CREATE INDEX DROAN_c ON Drosophila_ananassae(clust_id);
CREATE INDEX DROAN_s ON Drosophila_ananassae(species);
CREATE INDEX DROGR_c ON Drosophila_grimshawi(clust_id);
CREATE INDEX DROGR_s ON Drosophila_grimshawi(species);
CREATE INDEX DROME_c ON Drosophila_melanogaster(clust_id);
CREATE INDEX DROME_s ON Drosophila_melanogaster(species);
CREATE INDEX DROMO_c ON Drosophila_mojavensis(clust_id);
CREATE INDEX DROMO_s ON Drosophila_mojavensis(species);
CREATE INDEX DROPS_c ON Drosophila_pseudoobscura(clust_id);
CREATE INDEX DROPS_s ON Drosophila_pseudoobscura(species);
CREATE INDEX DROVI_c ON Drosophila_virilis(clust_id);
CREATE INDEX DROVI_s ON Drosophila_virilis(species);
CREATE INDEX DROWI_c ON Drosophila_willistoni(clust_id);
CREATE INDEX DROWI_s ON Drosophila_willistoni(species);
CREATE INDEX ENTHI_c ON Entamoeba_histolytica(clust_id);
CREATE INDEX ENTHI_s ON Entamoeba_histolytica(species);
CREATE INDEX EQUCA_c ON Equus_caballus(clust_id);
CREATE INDEX EQUCA_s ON Equus_caballus(species);
CREATE INDEX ESCCO_c ON Escherichia_coliK12(clust_id);
CREATE INDEX ESCCO_s ON Escherichia_coliK12(species);
CREATE INDEX FUSGR_c ON Fusarium_graminearum(clust_id);
CREATE INDEX FUSGR_s ON Fusarium_graminearum(species);
CREATE INDEX GALGA_c ON Gallus_gallus(clust_id);
CREATE INDEX GALGA_s ON Gallus_gallus(species);
CREATE INDEX GASAC_c ON Gasterosteus_aculeatus(clust_id);
CREATE INDEX GASAC_s ON Gasterosteus_aculeatus(species);
CREATE INDEX GIALA_c ON Giardia_lamblia(clust_id);
CREATE INDEX GIALA_s ON Giardia_lamblia(species);
CREATE INDEX HELRO_c ON Helobdella_robusta(clust_id);
CREATE INDEX HELRO_s ON Helobdella_robusta(species);
CREATE INDEX IXOSC_c ON Ixodes_scapularis(clust_id);
CREATE INDEX IXOSC_s ON Ixodes_scapularis(species);
CREATE INDEX KLULA_c ON Kluyveromyces_lactis(clust_id);
CREATE INDEX KLULA_s ON Kluyveromyces_lactis(species);
CREATE INDEX LEIMA_c ON Leishmania_major(clust_id);
CREATE INDEX LEIMA_s ON Leishmania_major(species);
CREATE INDEX LOTGI_c ON Lottia_gigantea(clust_id);
CREATE INDEX LOTGI_s ON Lottia_gigantea(species);
CREATE INDEX MACMU_c ON Macaca_mulatta(clust_id);
CREATE INDEX MACMU_s ON Macaca_mulatta(species);
CREATE INDEX MAGGR_c ON Magnaporthe_grisea(clust_id);
CREATE INDEX MAGGR_s ON Magnaporthe_grisea(species);
CREATE INDEX MONDO_c ON Monodelphis_domestica(clust_id);
CREATE INDEX MONDO_s ON Monodelphis_domestica(species);
CREATE INDEX MONBR_c ON Monosiga_brevicollis(clust_id);
CREATE INDEX MONBR_s ON Monosiga_brevicollis(species);
CREATE INDEX MUSMU_c ON Mus_musculus(clust_id);
CREATE INDEX MUSMU_s ON Mus_musculus(species);
CREATE INDEX NASVI_c ON Nasonia_vitripennis(clust_id);
CREATE INDEX NASVI_s ON Nasonia_vitripennis(species);
CREATE INDEX NEMVE_c ON Nematostella_vectensis(clust_id);
CREATE INDEX NEMVE_s ON Nematostella_vectensis(species);
CREATE INDEX NEUCR_c ON Neurospora_crassa(clust_id);
CREATE INDEX NEUCR_s ON Neurospora_crassa(species);
CREATE INDEX ORNAN_c ON Ornithorhynchus_anatinus(clust_id);
CREATE INDEX ORNAN_s ON Ornithorhynchus_anatinus(species);
CREATE INDEX ORYSA_c ON Oryza_sativa(clust_id);
CREATE INDEX ORYSA_s ON Oryza_sativa(species);
CREATE INDEX ORYLA_c ON Oryzias_latipes(clust_id);
CREATE INDEX ORYLA_s ON Oryzias_latipes(species);
CREATE INDEX OSTTA_c ON Ostreococcus_tauri(clust_id);
CREATE INDEX OSTTA_s ON Ostreococcus_tauri(species);
CREATE INDEX PANTR_c ON Pan_troglodytes(clust_id);
CREATE INDEX PANTR_s ON Pan_troglodytes(species);
CREATE INDEX PEDPA_c ON Pediculus_humanus(clust_id);
CREATE INDEX PEDPA_s ON Pediculus_humanus(species);
CREATE INDEX PHYPA_c ON Physcomitrella_patens(clust_id);
CREATE INDEX PHYPA_s ON Physcomitrella_patens(species);
CREATE INDEX PHYRA_c ON Phytophthora_ramorum(clust_id);
CREATE INDEX PHYRA_s ON Phytophthora_ramorum(species);
CREATE INDEX PHYSO_c ON Phytophthora_sojae(clust_id);
CREATE INDEX PHYSO_s ON Phytophthora_sojae(species);
CREATE INDEX PLAFA_c ON Plasmodium_falciparum(clust_id);
CREATE INDEX PLAFA_s ON Plasmodium_falciparum(species);
CREATE INDEX PLAVI_c ON Plasmodium_vivax(clust_id);
CREATE INDEX PLAVI_s ON Plasmodium_vivax(species);
CREATE INDEX PONPY_c ON Pongo_pygmaeus(clust_id);
CREATE INDEX PONPY_s ON Pongo_pygmaeus(species);
CREATE INDEX POPTR_c ON Populus_trichocarpa(clust_id);
CREATE INDEX POPTR_s ON Populus_trichocarpa(species);
CREATE INDEX PRIPA_c ON Pristionchus_pacificus(clust_id);
CREATE INDEX PRIPA_s ON Pristionchus_pacificus(species);
CREATE INDEX PUCGR_c ON Puccinia_graminis(clust_id);
CREATE INDEX PUCGR_s ON Puccinia_graminis(species);
CREATE INDEX RATNO_c ON Rattus_norvegicus(clust_id);
CREATE INDEX RATNO_s ON Rattus_norvegicus(species);
CREATE INDEX RHIOR_c ON Rhizopus_oryzae(clust_id);
CREATE INDEX RHIOR_s ON Rhizopus_oryzae(species);
CREATE INDEX SACCE_c ON Saccharomyces_cerevisiae(clust_id);
CREATE INDEX SACCE_s ON Saccharomyces_cerevisiae(species);
CREATE INDEX SCHMA_c ON Schistosoma_mansoni(clust_id);
CREATE INDEX SCHMA_s ON Schistosoma_mansoni(species);
CREATE INDEX SCHPO_c ON Schizosaccharomyces_pombe(clust_id);
CREATE INDEX SCHPO_s ON Schizosaccharomyces_pombe(species);
CREATE INDEX SCLSC_c ON Sclerotinia_sclerotiorum(clust_id);
CREATE INDEX SCLSC_s ON Sclerotinia_sclerotiorum(species);
CREATE INDEX SORBI_c ON Sorghum_bicolor(clust_id);
CREATE INDEX SORBI_s ON Sorghum_bicolor(species);
CREATE INDEX STANO_c ON Stagonospora_nodorum(clust_id);
CREATE INDEX STANO_s ON Stagonospora_nodorum(species);
CREATE INDEX STRPU_c ON Strongylocentrotus_purpuratus(clust_id);
CREATE INDEX STRPU_s ON Strongylocentrotus_purpuratus(species);
CREATE INDEX TAKRU_c ON Takifugu_rubripes(clust_id);
CREATE INDEX TAKRU_s ON Takifugu_rubripes(species);
CREATE INDEX TETTH_c ON Tetrahymena_thermophila(clust_id);
CREATE INDEX TETTH_s ON Tetrahymena_thermophila(species);
CREATE INDEX TETNI_c ON Tetraodon_nigroviridis(clust_id);
CREATE INDEX TETNI_s ON Tetraodon_nigroviridis(species);
CREATE INDEX THAPS_c ON Thalassiosira_pseudonana(clust_id);
CREATE INDEX THAPS_s ON Thalassiosira_pseudonana(species);
CREATE INDEX THEAN_c ON Theileria_annulata(clust_id);
CREATE INDEX THEAN_s ON Theileria_annulata(species);
CREATE INDEX THEPA_c ON Theileria_parva(clust_id);
CREATE INDEX THEPA_s ON Theileria_parva(species);
CREATE INDEX TRICA_c ON Tribolium_castaneum(clust_id);
CREATE INDEX TRICA_s ON Tribolium_castaneum(species);
CREATE INDEX TRIVA_c ON Trichomonas_vaginalis(clust_id);
CREATE INDEX TRIVA_s ON Trichomonas_vaginalis(species);
CREATE INDEX TRIAD_c ON Trichoplax_adhaerens(clust_id);
CREATE INDEX TRIAD_s ON Trichoplax_adhaerens(species);
CREATE INDEX TRYCR_c ON Trypanosoma_cruzi(clust_id);
CREATE INDEX TRYCR_s ON Trypanosoma_cruzi(species);
CREATE INDEX USTMA_c ON Ustilago_maydis(clust_id);
CREATE INDEX USTMA_s ON Ustilago_maydis(species);
CREATE INDEX XENTR_c ON Xenopus_tropicalis(clust_id);
CREATE INDEX XENTR_s ON Xenopus_tropicalis(species);
CREATE INDEX YARLI_c ON Yarrowia_lipolytica(clust_id);
CREATE INDEX YARLI_s ON Yarrowia_lipolytica(species);
