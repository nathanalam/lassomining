## CLI Version of Lasso Mining Algorithm

This is a simplified command line tool that simply reads FASTA information from Standard Input and writes to Standard Output the lasso peptide results.

## Prerequisites

Need to have MEME-SUITE installed.

## Examples

# Directly search from NCBI
esearch -db nucleotide -query "NZ_CP033730.1" | efetch -format fasta | python3 lassoMine.py > test.text