
# Section 1

# Section 1.1

# Load Biostrings if not already loaded
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings", update = FALSE)
  library(Biostrings)
} else {
  library(Biostrings)
}


fasta_file <- "/Users/jameswatson/Downloads/ncbi_dataset1/ncbi_dataset/data/GCF_000006765.1/cds_from_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")


# Section 1.2

# TRACTS COMPARISION

sequence_string <- paste(sequences, collapse = "")

# G tracts

count_motifs_G <- function(sequence, motif_G) {
  positions <- gregexpr(motif_G, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs_G <- c("GG", "GGG", "GGGG", "GGGGG")


motif_counts_G <- sapply(motifs_G, function(motif_G) {
  count_motifs_G(sequence_string, motif_G)
})

# Display counts
for (i in seq_along(motifs_G)) {
  cat(motifs_G[i], "motifs found:", motif_counts_G[i], "\n")
}

# C Tracts

count_motifs_C <- function(sequence, motif_C) {
  positions <- gregexpr(motif_C, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs_C <- c("CC", "CCC", "CCCC", "CCCCC")


motif_counts_C <- sapply(motifs_C, function(motif_C) {
  count_motifs_C(sequence_string, motif_C)
})

# Display counts
for (i in seq_along(motifs_C)) {
  cat(motifs_C[i], "motifs found:", motif_counts_C[i], "\n")
}

# A Tracts

count_motifs_A <- function(sequence, motif_A) {
  positions <- gregexpr(motif_A, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs_A <- c("AA", "AAA", "AAAA", "AAAAA")


motif_counts_A <- sapply(motifs_A, function(motif_A) {
  count_motifs_A(sequence_string, motif_A)
})

# Display counts
for (i in seq_along(motifs_A)) {
  cat(motifs_A[i], "motifs found:", motif_counts_A[i], "\n")
}

# T Tracts

count_motifs_T <- function(sequence, motif_T) {
  positions <- gregexpr(motif_T, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs_T <- c("TT", "TTT", "TTTT", "TTTTT")


motif_counts_T <- sapply(motifs_T, function(motif_T) {
  count_motifs_T(sequence_string, motif_T)
})

# Display counts
for (i in seq_along(motifs_T)) {
  cat(motifs_T[i], "motifs found:", motif_counts_T[i], "\n")
}

# Plotting and comparing these tracts 

library(ggplot2)

# Create a data frame 
motif_data <- data.frame(
  motifs = rep(c("XX", "XXX", "XXXX", "XXXXX"), 4),
  bases = rep(c("A", "G", "C", "T"), each = 4),
  counts = c(motif_counts_A, motif_counts_G, motif_counts_C, motif_counts_T)
)


motif_data <- motif_data[order(match(motif_data$motifs, c("XX", "XXX", "XXXX", "XXXXX")), motif_data$bases), ]

# Plotting
ggplot(motif_data, aes(x = motifs, y = counts, fill = bases)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = counts), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(x = "Motifs", y = "Frequency", title = "Frequency of Motifs by Base") +
  scale_fill_manual(values = c("A" = "blueviolet", "G" = "aquamarine", "C" = "maroon2", "T" = "steelblue3")) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

#-----------------------------------------------------------------

# Section 1.3

# Homing in on triplet G tracts (AGGGA, AGGGT, AGGGC | CGGGC, CGGGA, CGGGT | TGGGT, TGGGC, TGGGA)



count_motifs <- function(sequence, motif) {
  positions <- gregexpr(motif, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs <- c("AGGGTA", "AGGGTT", "AGGGTC", "AGGGTG", "CGGGTC", "CGGGTA", "CGGGTT", "CGGGTG", "TGGGTT", "TGGGTC", "TGGGTA", "TGGGTG")


motif_counts <- sapply(motifs, function(motif) {
  count_motifs(sequence_string, motif)
})


motif_data <- data.frame(
  motifs = motifs,
  counts = motif_counts,
  bases = substr(motifs, 1, 1)  
)

# Plotting
library(ggplot2)

ggplot(motif_data, aes(x = motifs, y = counts)) +
  geom_bar(stat = "identity", fill = "blueviolet") +
  geom_text(aes(label = counts), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(x = "Motifs", y = "Frequency", title = "Frequency of Triplet G Tract Motifs") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")



# Homing in on quad G tracts (AGGGGA, AGGGGT, AGGGGC | CGGGGC, CGGGGA, CGGGGT | TGGGGT, TGGGGC, TGGGGA)



count_motifs <- function(sequence, motif) {
  positions <- gregexpr(motif, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs <- c("AGGGGTA", "AGGGGTT", "AGGGGTC", "AGGGGTG", "CGGGGTC", "CGGGGTA", "CGGGGTT", "CGGGGTG", "TGGGGTT", "TGGGGTC", "TGGGGTA", "TGGGGTG")


motif_counts <- sapply(motifs, function(motif) {
  count_motifs(sequence_string, motif)
})


motif_data <- data.frame(
  motifs = motifs,
  counts = motif_counts,
  bases = substr(motifs, 1, 1)  
)

# Plotting
library(ggplot2)

ggplot(motif_data, aes(x = motifs, y = counts)) +
  geom_bar(stat = "identity", fill = "blueviolet") +
  geom_text(aes(label = counts), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(x = "Motifs", y = "Frequency", title = "Frequency of Quad G Tract Motifs") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")



# Homing in on Pentet G tracts (AGGGGGA, AGGGGGT, AGGGGGC | CGGGGGC, CGGGGGA, CGGGGGT | TGGGGGT, TGGGGGC, TGGGGGA)



count_motifs <- function(sequence, motif) {
  positions <- gregexpr(motif, sequence)
  count <- sum(unlist(positions) != -1)
  return(count)
}


motifs <- c("AGGGGGTA", "AGGGGGTT", "AGGGGGTC", "AGGGGGTG", "CGGGGGTC", "CGGGGGTA", "CGGGGGTT", "CGGGGGTG", "TGGGGGTT", "TGGGGGTC", "TGGGGGTA", "TGGGGGTG")


motif_counts <- sapply(motifs, function(motif) {
  count_motifs(sequence_string, motif)
})


motif_data <- data.frame(
  motifs = motifs,
  counts = motif_counts,
  bases = substr(motifs, 1, 1)  
)

# Plotting
library(ggplot2)

ggplot(motif_data, aes(x = motifs, y = counts)) +
  geom_bar(stat = "identity", fill = "blueviolet") +
  geom_text(aes(label = counts), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(x = "Motifs", y = "Frequency", title = "Frequency of Pentet G Tract Motifs") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")



# Section 1.4

# Percentage nucleotide counts 

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")


total_nucleotides <- nchar(gsub("[^ACGT]", "", sequence_string))


a_count <- nchar(gsub("[^A]", "", sequence_string))
c_count <- nchar(gsub("[^C]", "", sequence_string))
g_count <- nchar(gsub("[^G]", "", sequence_string))
t_count <- nchar(gsub("[^T]", "", sequence_string))


a_percentage <- (a_count / total_nucleotides) * 100
c_percentage <- (c_count / total_nucleotides) * 100
g_percentage <- (g_count / total_nucleotides) * 100
t_percentage <- (t_count / total_nucleotides) * 100


cat("Percentage of A:", a_percentage, "%\n")
cat("Percentage of C:", c_percentage, "%\n")
cat("Percentage of G:", g_percentage, "%\n")
cat("Percentage of T:", t_percentage, "%\n")






# Section 2


# Section 2.1

# Positions and locations of AGGGGGTC (pentet) motifs 


library(GenomicRanges)
library(rtracklayer)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

motif <- "AGGGGGTC"

matches <- matchPattern(motif, DNAString(sequence_string))

positions <- start(matches)

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


for (X in positions) {
  
  position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
  overlaps <- findOverlaps(position_gr, gff_data)
  
  if (length(overlaps) > 0) {
    overlapping_features <- gff_data[subjectHits(overlaps)]
    
    
    valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
    
    if (length(valid_features) > 0) {
      locus_tags <- mcols(valid_features)$locus_tag
      symbols <- mcols(valid_features)$gene
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      unique_symbols <- unique(na.omit(symbols))
      
      final_locus_tag <- if(length(unique_locus_tags) > 0) unique_locus_tags[1] else "Unknown"
      final_symbol <- if(length(unique_symbols) > 0) unique_symbols[1] else "Unknown"
    } else {
      final_locus_tag <- "Unknown"
      final_symbol <- "Unknown"
    }
    
    cat("Position:", X, "- Locus Tag:", final_locus_tag, "Symbol:", final_symbol, "\n")
  } else {
    cat("Position:", X, "- No features overlap with the specified position.\n")
  }
}



# Section 2.2

# Amino acid mutant sequence encoded via leading strand genes 

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)


codon_to_amino_acid <- list(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "X", TAG = "X",
  TGT = "C", TGC = "C", TGA = "X", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)


translateSequence <- function(nucleotide_seq) {
  
  codons <- strsplit(nucleotide_seq, split = "")[[1]]
  codons <- paste(codons[seq(1, length(codons), 3)], 
                  codons[seq(2, length(codons), 3)], 
                  codons[seq(3, length(codons), 3)], 
                  sep = "")
  
  amino_acids <- sapply(codons, function(codon) if(!is.null(codon_to_amino_acid[[codon]])) codon_to_amino_acid[[codon]] else "-")
  
  return(paste(amino_acids, collapse = ""))
}


fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


locus_tags <- c("PA0138", "PA0222", "PA0828", "PA0970", "PA1048", "PA1180", "PA1205", "PA1742", "PA1974", "PA2485", "PA2638", "PA3837", "PA3965", "PA4321", "PA4675")


genes_of_interest <- gff_data[gff_data$type == "gene" & gff_data$locus_tag %in% locus_tags, ]


for (i in seq_along(genes_of_interest$locus_tag)) {
  start <- start(genes_of_interest)[i]
  end <- end(genes_of_interest)[i]
  sequence <- substr(sequence_string, start, end)
  
  
  mutated_sequence <- gsub("AGGGGGTC", "AGGGGGGC", sequence)
  
  
  original_aa_sequence <- translateSequence(sequence)
  mutated_aa_sequence <- translateSequence(mutated_sequence)
  
  
  amino_acid_mutations <- findAminoAcidMutations(original_aa_sequence, mutated_aa_sequence)
  
  cat(sprintf("Locus Tag: %s\nOriginal Sequence: %s\nMutated Sequence: %s\nOriginal Amino Acid Sequence: %s\nMutated Amino Acid Sequence: %s\nAmino Acid Mutations: %s\n\n",
              genes_of_interest$locus_tag[i], sequence, mutated_sequence, original_aa_sequence, mutated_aa_sequence, paste(sapply(amino_acid_mutations, function(x) paste(x, collapse=", ")), collapse=", ")))
}


# Section 2.3 

# Amino acid mutant sequence encoded via lagging strand genes 

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)


codon_to_amino_acid <- list(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "X", TAG = "X",
  TGT = "C", TGC = "C", TGA = "X", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)


translateSequence <- function(nucleotide_seq) {
  codons <- strsplit(nucleotide_seq, split = "")[[1]]
  codons <- paste(codons[seq(1, length(codons), 3)], 
                  codons[seq(2, length(codons), 3)], 
                  codons[seq(3, length(codons), 3)], 
                  sep = "")
  amino_acids <- sapply(codons, function(codon) if(!is.null(codon_to_amino_acid[[codon]])) codon_to_amino_acid[[codon]] else "-")
  return(paste(amino_acids, collapse = ""))
}


reverseComplementSequence <- function(nucleotide_seq) {
  reversed_seq <- rev(strsplit(nucleotide_seq, "")[[1]])
  complemented_seq <- sapply(reversed_seq, function(nucleotide) {
    if(nucleotide == "A") return("T")
    else if(nucleotide == "T") return("A")
    else if(nucleotide == "C") return("G")
    else if(nucleotide == "G") return("C")
    else return(nucleotide)
  })
  return(paste(complemented_seq, collapse = ""))
}


fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


locus_tags <- c("PA0293", "PA1044", "PA1107", "PA1779", "PA2906", "PA3109", "PA3596", "PA3769", "PA3770", "PA3804", "PA3943", "PA3950", "PA4812", "PA5292")


genes_of_interest <- gff_data[gff_data$type == "gene" & gff_data$locus_tag %in% locus_tags, ]


findAminoAcidMutations <- function(original_aa_seq, mutated_aa_seq) {
  mutations <- list()
  for (i in 1:nchar(original_aa_seq)) {
    original_aa <- substr(original_aa_seq, i, i)
    mutated_aa <- substr(mutated_aa_seq, i, i)
    if (original_aa != mutated_aa) {
      mutations[[length(mutations) + 1]] <- sprintf("Position %d: %s -> %s", i, original_aa, mutated_aa)
    }
  }
  return(mutations)
}


for (i in seq_along(genes_of_interest$locus_tag)) {
  start <- start(genes_of_interest)[i]
  end <- end(genes_of_interest)[i]
  sequence <- substr(sequence_string, start, end)
  
  
  original_rev_compl_sequence <- reverseComplementSequence(sequence)
  
  
  original_aa_sequence <- translateSequence(original_rev_compl_sequence)
  
  
  mutated_sequence <- gsub("AGGGGGTC", "AGGGGGGC", sequence)
  
  
  mutated_rev_compl_sequence <- reverseComplementSequence(mutated_sequence)
  
  
  mutated_aa_sequence <- translateSequence(mutated_rev_compl_sequence)
  
  
  amino_acid_mutations <- findAminoAcidMutations(original_aa_sequence, mutated_aa_sequence)
  
  cat(sprintf("Locus Tag: %s\nOriginal Sequence: %s\nMutated Sequence: %s\nOriginal Reversed & Complemented Sequence: %s\nMutated Reversed & Complemented Sequence: %s\nOriginal Amino Acid Sequence: %s\nMutated Amino Acid Sequence: %s\nAmino Acid Mutations: %s\n\n",
              genes_of_interest$locus_tag[i], sequence, mutated_sequence, original_rev_compl_sequence, mutated_rev_compl_sequence, original_aa_sequence, mutated_aa_sequence, paste(sapply(amino_acid_mutations, function(x) paste(x, collapse=", ")), collapse=", ")))
}



# Section 2.4

#Bootstrap test to determine if secondary structures are enriched or depleted

array <- matrix(nrow=100000, ncol = 2)

#nrow == how many times you want to repeat the bootstrap. Start with 100k and see how long that takes, then scale up. Aim for 500k-1mil if you can.

#ncol == how many outcomes there are. This has to be left at two for this pipeline to work.

weight <- c(0.667,0.333)

#weight == how likely each of the two outcomes are to happen. If you're interested in say one mutation in a set of 20, then you'd assign a likelihood to the one mutation and then sum the other 19 for the other column. These should always sum to 1.

for(i in 1:100000){ #The number after 1:n should match nrow
  
  dis <- sample(1:2, 30, replace=TRUE, prob=weight)
  
  #The second number listed here describes how many 'draws' you're doing per iteration. For example,if you did an experiment where you sequenced 21 samples and saw x mutation 4 times, you'd set this number to 21.
  
  distab <- as.data.frame(table(dis))
  
  array[i, 1] = distab[2, 2]
  
}

array[is.na(array)] <- 0 #replaces na values with 0

array <- array[order(array[,2]),] #Organises array so that it's in ascending order relative to column 3


summary <- as.data.frame(table(array))
print(summary)


frequencies <- c(8, 60, 275, 921, 2329, 4846, 8372,
                 11927, 14476, 15135, 13868, 11054, 7763,
                 4555, 2463, 1204, 473, 185, 71, 9, 4, 1)


total_samples <- sum(frequencies)


p_value_depletion <- sum(frequencies[1:4]) / total_samples


p_value_enrichment <- sum(frequencies[18:length(frequencies)]) / total_samples

cat("P-value for depletion (4 or fewer out of 30):", p_value_depletion, "\n")
cat("P-value for enrichment (19 or more out of 30):", p_value_enrichment, "\n")



# Section 3

                 
# Are G tracks enriched in particular gene classes? 

# Section 3.1

## Quad G tracks 

### Two-Component systems 

setwd("/Users/jameswatson/Downloads")
two_component_genes <- read.csv("features.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, two_component_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within two-component system genes:", length(matching_genes), "\n")



### Transcriptional regulators 

setwd("/Users/jameswatson/Downloads")
transcriptional_regulator_genes <- read.csv("features (1).csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, transcriptional_regulator_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within transcriptional regulator genes:", length(matching_genes), "\n")


### Motility and attachment 

setwd("/Users/jameswatson/Downloads")
motility_attachment_genes <- read.csv("features (2).csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, motility_attachment_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within motility and attachment genes:", length(matching_genes), "\n")



### Cell wall / LPS / capsule 

setwd("/Users/jameswatson/Downloads")
cell_wall_LPS_capsule_genes <- read.csv("features (3).csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, cell_wall_LPS_capsule_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within cell wall / LPS / capsule genes:", length(matching_genes), "\n")



### Transport of small molecules 

setwd("/Users/jameswatson/Downloads")
transport_of_small_molecule_genes <- read.csv("features (4).csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, transport_of_small_molecule_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within transport of small molecule genes:", length(matching_genes), "\n")


### AMR

setwd("/Users/jameswatson/Downloads")
AMR_genes <- read.csv("AMR_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, AMR_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within AMR genes:", length(matching_genes), "\n")


### Membrane proteins 

setwd("/Users/jameswatson/Downloads")
Membrane_proteins_genes <- read.csv("Membrane_proteins_genes.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)


motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA", "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA", "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


all_genes_with_motifs <- character(0)

for (motif in motifs) {
  matches <- matchPattern(motif, DNAString(sequence_string))
  positions <- start(matches)
  
  for (X in positions) {
    position_gr <- GRanges(seqnames = "NC_002516.2", ranges = IRanges(start = X, end = X + nchar(motif) - 1))
    overlaps <- findOverlaps(position_gr, gff_data)
    
    if (length(overlaps) > 0) {
      overlapping_features <- gff_data[subjectHits(overlaps)]
      valid_features <- overlapping_features[which(start(overlapping_features) <= X & end(overlapping_features) >= X + nchar(motif) - 1)]
      locus_tags <- mcols(valid_features)$locus_tag
      
      unique_locus_tags <- unique(na.omit(locus_tags))
      all_genes_with_motifs <- unique(c(all_genes_with_motifs, unique_locus_tags))
    }
  }
}


matching_genes <- intersect(all_genes_with_motifs, Membrane_proteins_genes$Locus.Tag)


cat("Number of G4T motif containing genes found within Membrane protein genes:", length(matching_genes), "\n")


#### Plotting results (run previous codes first)

library(ggplot2)


n_two_component_system <- length(intersect(all_genes_with_motifs, two_component_genes$Locus.Tag))
n_transcriptional_regulators <- length(intersect(all_genes_with_motifs, transcriptional_regulator_genes$Locus.Tag))
n_motility_attachment <- length(intersect(all_genes_with_motifs, motility_attachment_genes$Locus.Tag))
n_cell_wall_LPS_capsule <- length(intersect(all_genes_with_motifs, cell_wall_LPS_capsule_genes$Locus.Tag))
n_transport_of_small_molecule <- length(intersect(all_genes_with_motifs, transport_of_small_molecule_genes$Locus.Tag))
n_AMR <- length(intersect(all_genes_with_motifs, AMR_genes$Locus.Tag))
n_membrane_proteins <- length(intersect(all_genes_with_motifs, Membrane_proteins_genes$Locus.Tag))


data <- data.frame(
  Category = c("Two-component system", "Transcriptional Regulators", "Motility & Attachment", "Cell Wall/LPS/Capsule", "Transport of Small Molecules", "Antibiotic Resistance & Susceptibility", "Membrane Proteins"),
  Count = signif(c(n_two_component_system / 123 * 100, 
                   n_transcriptional_regulators / 492 * 100, 
                   n_motility_attachment / 140 * 100, 
                   n_cell_wall_LPS_capsule / 193 * 100, 
                   n_transport_of_small_molecule / 609 * 100, 
                   n_AMR / 74 * 100, 
                   n_membrane_proteins / 675 * 100), digits = 3))


# Plot
ggplot(data, aes(x = Category, y = Count)) +
  geom_bar(stat = "identity", fill = "grey18", show.legend = FALSE) +  
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.1, vjust = 0.5) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "Percentage of G4T Motifs Across Gene Classes",
       x = "Gene Classes",
       y = "Probability of a gene containing a G4T Motif") +
  theme(text = element_text(size = 12),
        plot.title = element_text(face = "bold"))



# Section 3.2

# Stats - Now we want to write a pipeline that generates a dummy dataset to see how likely G4T motifs are found in a set of genes by chance

# Extended with loop counter 

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

fasta_file <- "/Users/jameswatson/Downloads/2.PA01/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
sequences <- readDNAStringSet(fasta_file)
sequence_string <- paste(sequences, collapse = "")

gff_path <- "/Users/jameswatson/Downloads/PA01_genome_annotation/ncbi_dataset/data/GCF_000006765.1/genomic.gff" 
gff_data <- import(gff_path)

motifs <- c("AGGGGTC", "AGGGGTG", "AGGGGTT", "AGGGGTA",
            "TGGGGTC", "TGGGGTG", "TGGGGTT", "TGGGGTA",
            "CGGGGTC", "CGGGGTG", "CGGGGTT", "CGGGGTA")


total_motifs_counts <- numeric(10000)


total_iterations <- 10000


for (j in 1:total_iterations) {
  
  genes_of_interest <- subset(gff_data, type == "gene")
  selected_indices <- sample(length(genes_of_interest), 311)
  selected_genes <- genes_of_interest[selected_indices]
  
  
  extracted_sequences <- sapply(seq_along(selected_genes), function(i) {
    start <- start(selected_genes)[i]
    end <- end(selected_genes)[i]
    as.character(subseq(sequence_string, start, end))
  })
  concatenated_sequence <- paste(extracted_sequences, collapse = "")
  
  
  motif_counts <- sapply(motifs, function(motif) {
    sum(gregexpr(motif, concatenated_sequence, fixed = TRUE)[[1]] > 0)
  })
  total_motifs_counts[j] <- sum(motif_counts)
  
  
  if (j %% 100 == 0) {
    cat("Completed", j, "out of", total_iterations, "iterations.\n")
  }
}


average_motifs_found <- mean(total_motifs_counts)
cat("Average motifs found over 10,000 iterations:", average_motifs_found, "\n")


# Histogram plot with normality using just percentages 

percentages_per_gene <- (total_motifs_counts / 311) * 100


percentages_range <- range(percentages_per_gene)
print(percentages_range)

breaks_vector <- seq(30, 85, by = 5)


hist(percentages_per_gene, breaks = breaks_vector, probability = TRUE, 
     main = "Null Histogram Distribution", 
     xlab = "Percentage of G4T Motifs per Gene", ylab = "Density", col = "gray", border = "black", 
     ylim = c(0, max(density(percentages_per_gene)$y) * 1))


dens <- density(percentages_per_gene)
normal_dens <- curve(dnorm(x, mean = mean(percentages_per_gene), sd = sd(percentages_per_gene)), col = "blue", lwd = 2, add = TRUE)


abline(v = 57.2010932, col = "red", lwd = 2, lty = 2)
abline(v = 48.3, col = "darkmagenta", lwd = 2, lty = 2)
abline(v = 46.6, col = "royalblue", lwd = 2, lty = 2)
abline(v = 38.6, col = "turquoise1", lwd = 2, lty = 2)
abline(v = 31, col = "hotpink1", lwd = 2, lty = 2)
abline(v = 39.8, col = "slateblue2", lwd = 2, lty = 2)
abline(v = 51.4, col = "seagreen2", lwd = 2, lty = 2)
abline(v = 42.2, col = "orangered", lwd = 2, lty = 2)

legend("right", 
       legend = c("Normal Fit", "Null Average (57.2%)", "Transport of Small Molecules (48.3%)", 
                  "Cell Wall/LPS/Capsule (46.6%)", "Motility & Attachment (38.6%)", 
                  "Transcriptional Regulators (31.0%)", "Two-Component Systems (39.8%)", 
                  "Antibiotic Resistance & Susceptibility (51.4%)", "Membrane Proteins (42.2%)"),
       col = c("blue", "red", "darkmagenta", "royalblue", "turquoise1", 
               "hotpink1", "slateblue2", "seagreen2", "orangered"),
       lty = c(1, 2, 2, 2, 2, 2, 2, 2, 2),
       cex = 0.75,  
       xpd = TRUE,  
       inset = c(-0.25, 0),  
       bty = "n") 


# Section 3.3

# Z-test

percentages_per_gene <- (total_motifs_counts / 311) 

mean_percentage <- mean(percentages_per_gene)

std_dev <- sd(percentages_per_gene)

observed_value <- 0.422

z_score <- (observed_value - mean_percentage) / std_dev

p_value <- pnorm(z_score)

cat("Z-Score:", z_score, "\n")
cat("P-Value:", p_value, "\n")












