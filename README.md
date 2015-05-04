# HTS_codon_analyzer
 Reads, filters and analyzes high-throughput sequencing (HTS) data 
 
 
 Created by AGNES Toth-Petroczy, 31.05.2011
 Last modified 23.11.2013

## Input: alignment of the Illumina sequencing reads, WT sequence (fasta), gene_position_start, gene_position_end
 (aligned by http://www.novocraft.com/products/novoalign/)
## Output: 
 *.allmut (all mutations per codons per position)
 *.freq (frequencies of mutations (number of occurences/number of reads))
 *.ns1 (single numcleotide changes)
 *.ns2 (double nucleotide changes)
 *.ns3 (triple nucleotide changes),
 *.syn (synonymous changes)
 *.indels (insertions and deletions)
 *.codoposinreads (position within the reads (3-37 out of 40 bp long reads) where the mutation occurs)

##Example run
 ./HTS_filtered_codon_mutations.pl novo.align WTgene 31 1036
