This directory contains version 0.9 of the full set of schemas used in the
sqlite-based annotation data packages.

Version 0.9 is our target for the BioC 2.1 release i.e. it is the version
that we plan to use in the data packages that will be released with
Bioconductor 2.1 (AnnotationDbi 1.0.0).

Version 0.9 has been successfully tested (i.e. imported) with SQLite
(3.4.1), MySQL+InnoDB (5.0.26) and PostgreSQL (8.1.9) on a 64-bit openSUSE
10.2 system. It has not been tested on Oracle yet.

All the *.sqlite files using one of the 0.9 schemas must set DBSCHEMAVERSION
to 0.9 in their "metadata" table.

See the DataTypes.txt file for all the data types used across the 0.9 schemas.

