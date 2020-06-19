This directory contains version 2.0 of the full set of schemas used in the
sqlite-based annotation data packages.

Version 2.0 is our target for the BioC 2.5 release.

All the *.sqlite files using one of the 2.0 schemas must set DBSCHEMAVERSION
to 2.0 in their "metadata" table.

See the DataTypes.txt file for all the data types used across the 2.0 schemas.


SUMMARY OF CHANGES SINCE VERSION 1.0
------------------------------------

Note that some of these changes are disruptive (i.e. they potentially break the
SQL queries written for version 1.0).

All probe-based schemas:

  o dropped all tables except the following: genbank, map_counts, map_metadata, metadata, probes
  o _id in the probes table is now replaced by a gene_id which is a foreign key to a matching org package
  o accessions table has been added to hold accessions associated with probes by the manufacturer

