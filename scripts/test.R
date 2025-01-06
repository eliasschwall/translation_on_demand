# Load necessary libraries
library(biomaRt)
library(dplyr)
library(tidyr)

# Step 1: Load the processed proteomic and transcriptomic data
proteomic_data_cleaned <- read.csv("/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/proteomic_data_cleaned.csv", check.names = FALSE)
transcriptomic_data_cleaned <- read.csv("/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/transcriptomic_data_cleaned.csv", check.names = FALSE)

# Step 2: Get gene exon lengths using Ensembl BioMart
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get exon data for all genes: exon starts, exon ends, and gene IDs
exon_data <- getBM(attributes = c('ensembl_gene_id', 'exon_chrom_start', 'exon_chrom_end'), mart = ensembl)

# Calculate unique exon lengths for each gene
gene_exon_lengths <- exon_data %>%
  rowwise() %>%
  mutate(exon_positions = list(seq(exon_chrom_start, exon_chrom_end))) %>%
  ungroup() %>%
  dplyr::select(ensembl_gene_id, exon_positions) %>%
  unnest(cols = c(exon_positions)) %>%
  distinct(ensembl_gene_id, exon_positions) %>%
  group_by(ensembl_gene_id) %>%
  summarise(Exonic_Length = n())

# Save the gene exon lengths to a CSV file
write.csv(gene_exon_lengths, "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/reference/mouse_gene_exon_lengths.csv", row.names = FALSE)


gene_exon_lengths <- read.csv("/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/reference/mouse_gene_exon_lengths.csv")
# Step 3: TPM Normalization
transcriptomic_data_cleaned <- merge(transcriptomic_data_cleaned, gene_exon_lengths, by = "ensembl_gene_id")
transcriptomic_data_cleaned$Exonic_Length <- transcriptomic_data_cleaned$Exonic_Length / 1000
row.names(transcriptomic_data_cleaned) <- transcriptomic_data_cleaned$ensembl_gene_id

# Normalize transcript counts by exon length (RPK)
transcriptomic_data_cleaned[, 3:(ncol(transcriptomic_data_cleaned) - 1)] <- transcriptomic_data_cleaned[, 3:(ncol(transcriptomic_data_cleaned) - 1)] / transcriptomic_data_cleaned$Exonic_Length
RPK <- transcriptomic_data_cleaned[, 3:(ncol(transcriptomic_data_cleaned) - 1)]

# TPM scaling factors
TPM_scaling_factors <- colSums(RPK) / 1e6

# TPM calculation
TPM <- sweep(RPK, 2, TPM_scaling_factors, "/")

# Log2 transformation
TPM_log2 <- log2(TPM + 1)

# Save the TPM and log2 TPM data
write.csv(RPK, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/transcriptomic_data_RPK.csv")
write.csv(TPM, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/transcriptomic_data_TPM.csv")
write.csv(TPM_log2, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/transcriptomic_data_TPM_log2.csv")

# Append GeneID and Ensembl ID columns
TPM_log2$GeneID <- transcriptomic_data_cleaned$GeneID
TPM_log2$ensembl_gene_id <- transcriptomic_data_cleaned$ensembl_gene_id

# Step 4: Merging proteomic and transcriptomic data
proteome_and_transcriptome_TPM_log2 <- merge(proteomic_data_cleaned, TPM_log2, by = "GeneID") %>%
  relocate(ensembl_gene_id, .after = GeneID) %>%
  setNames(c("GeneID", "EnsemblID", "UniprotID", "0m_protein", "1h_protein", "6h_protein", "12h_protein", "24h_protein", "36h_protein", "48h_protein", "72h_protein",
             "0m_transcript", "1h_transcript", "6h_transcript", "12h_transcript", "24h_transcript", "36h_transcript", "48h_transcript", "72h_transcript"))

# Save the merged data
write.csv(proteome_and_transcriptome_TPM_log2, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/proteome_and_transcriptome_TPM_log2.csv")

# Step 5: Fold Change (FC) Calculation Function
calculating_log2FC_between_time_points <- function(count_matrix) {
  # Extract the time points from the column names
  time_points <- colnames(count_matrix)
  
  # Initialize an empty dataframe with the same number of rows as count_matrix
  FC_df <- data.frame(matrix(ncol = 0, nrow = nrow(count_matrix)))
  
  # Calculate log2 fold changes (difference between consecutive log2-transformed values)
  for(i in 2:ncol(count_matrix)) {
    # Calculate log2 fold change as the difference between time points (since data is already log2-transformed)
    FC_column <- count_matrix[, i] - count_matrix[, i-1]
    
    # Add the log2 fold-change column to FC_df with the appropriate name
    FC_df[[paste0("log2FC_", time_points[i-1], "_to_", time_points[i])]] <- FC_column
  }
  
  # Set row names of FC_df to match the original count_matrix
  rownames(FC_df) <- rownames(count_matrix)
  
  return(FC_df)
}

# Step 6: Calculate log2 fold changes for transcriptomic and proteomic data
TPM_log2_FC <- calculating_log2FC_between_time_points(TPM_log2[, 1:8])
TPM_log2_FC$GeneID <- TPM_log2$GeneID 

proteomic_data_FC <- calculating_FC_between_time_points(proteomic_data_cleaned[, 3:10])
proteomic_data_FC$GeneID <- proteomic_data_cleaned$GeneID

# Save fold change data
write.csv(TPM_log2_FC, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/TPM_log2_FC.csv")
write.csv(proteomic_data_FC, file = "/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/processed/proteomic_data_FC.csv")