#!/bin/bash
#
# conda install -c bioconda entrez-direct --yes
#
# H1
efetch -id CM066638.1 -db nucleotide -format fasta > H1_C1.fa
efetch -id CM066639.1 -db nucleotide -format fasta > H1_C2.fa
efetch -id CM066640.1 -db nucleotide -format fasta > H1_C3.fa
efetch -id CM066641.1 -db nucleotide -format fasta > H1_C4.fa
efetch -id CM066642.1 -db nucleotide -format fasta > H1_C5.fa
efetch -id CM066643.1 -db nucleotide -format fasta > H1_C6.fa
efetch -id CM066644.1 -db nucleotide -format fasta > H1_C7.fa
efetch -id CM066645.1 -db nucleotide -format fasta > H1_C8.fa
efetch -id CM066646.1 -db nucleotide -format fasta > H1_C9.fa
efetch -id CM066647.1 -db nucleotide -format fasta > H1_C10.fa
efetch -id CM066648.1 -db nucleotide -format fasta > H1_C11.fa
efetch -id CM066649.1 -db nucleotide -format fasta > H1_C12.fa
efetch -id CM066650.1 -db nucleotide -format fasta > H1_C13.fa
efetch -id CM066651.1 -db nucleotide -format fasta > H1_C14.fa
efetch -id CM066652.1 -db nucleotide -format fasta > H1_C15.fa
efetch -id CM066653.1 -db nucleotide -format fasta > H1_C16.fa
efetch -id CM066654.1 -db nucleotide -format fasta > H1_C17.fa
efetch -id CM066655.1 -db nucleotide -format fasta > H1_C18.fa
# H2
efetch -id CM066656.1 -db nucleotide -format fasta > H2_C1.fa
efetch -id CM066657.1 -db nucleotide -format fasta > H2_C2.fa
efetch -id CM066658.1 -db nucleotide -format fasta > H2_C3.fa
efetch -id CM066659.1 -db nucleotide -format fasta > H2_C4.fa
efetch -id CM066660.1 -db nucleotide -format fasta > H2_C5.fa
efetch -id CM066661.1 -db nucleotide -format fasta > H2_C6.fa
efetch -id CM066662.1 -db nucleotide -format fasta > H2_C7.fa
efetch -id CM066663.1 -db nucleotide -format fasta > H2_C8.fa
efetch -id CM066664.1 -db nucleotide -format fasta > H2_C9.fa
efetch -id CM066665.1 -db nucleotide -format fasta > H2_C10.fa
efetch -id CM066666.1 -db nucleotide -format fasta > H2_C11.fa
efetch -id CM066667.1 -db nucleotide -format fasta > H2_C12.fa
efetch -id CM066668.1 -db nucleotide -format fasta > H2_C13.fa
efetch -id CM066669.1 -db nucleotide -format fasta > H2_C14.fa
efetch -id CM066670.1 -db nucleotide -format fasta > H2_C15.fa
efetch -id CM066671.1 -db nucleotide -format fasta > H2_C16.fa
efetch -id CM066672.1 -db nucleotide -format fasta > H2_C17.fa
efetch -id CM066673.1 -db nucleotide -format fasta > H2_C18.fa
#

cat H1_C*.fa > H1.fa
cat H2_C*.fa > H2.fa
rm -f H1_C*.fa H2_C*.fa
cat H1.fa H2.fa > TME117.fa
rm -f H1.fa H2.fa