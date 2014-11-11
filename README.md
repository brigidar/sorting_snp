sorting_snp
after snp discovery a costom script verifies for the positions to be "real" snps. The table has to be cleaned for positions that were wrongly called, such as misalignments in repeated regions, or positions that are not present in all query genomes. Table is read as a dataframe and only the bases are further processed to create fasta files for phylogeny.

===========
