# -----------------------------------------------------------------------------
#
# Compare ortholog lengths to percent identities
#
# Jason Jiang - Created: Apr/13/2023
#               Last edited: August/04/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(seqinr))

################################################################################

## Global variables

# Parameters for Needleman-Wunsch alignment
GAP_OPEN = 10.0
GAP_EXTEND = 0.5
ALIGNMENT_TYPE = 'global'
EXCLUDED_SP = c("N_bomb", "P_neur", "D_muel", "H_magn", "C_dike", "H_tvae", "N_apis", "D_roes")

################################################################################

main <- function() {
  ref_essential_genes <- read_lines('../../data/essential_yeast_genes.txt')
  clades <- get_clades_hash(read_csv('../../data/species_clades.csv', show_col_types = FALSE))
  
  orthologs <- read_csv('../../results/single_copy_orthogroups.csv') %>%
    filter(!(species %in% EXCLUDED_SP)) %>%
    mutate(essential_in_ref = ref_ortholog %in% ref_essential_genes) %>%
    rowwise() %>%
    mutate(clade = clades[[species]]) %>%
    ungroup() %>%
    mutate(species_is_outgroup = ifelse(clade == 'Outgroup', TRUE, FALSE))
  
  # create hashmap of ortholog sequences in orthogroups
  orthogroups_of_interest <- unique(orthologs$orthogroup)
  orthogroup_seq_hash <- new.env()
  
  for (og in orthogroups_of_interest) {
    fasta <- read.fasta(
      str_c('../../results/OrthoFinder/Results_OrthoFinder/Orthogroup_Sequences/',
            og, '.fa'),
      seqtype = 'AA',
      as.string = TRUE
    )
    
    orthogroup_seqs <- new.env()
    for (ortholog in getName(fasta)) {
      # replace any unconventional amino acids (ex: U) with their conventional
      # equivalents (ex: U -> C)
      orthogroup_seqs[[ortholog]] <- str_c(
        recode(str_split(fasta[[ortholog]][1], '')[[1]],
               !!!c(X = 'A', U = 'C', O = 'K', B = 'D', Z = 'E', J = 'L')),
        collapse = '')
    }
    
    orthogroup_seq_hash[[og]] <- orthogroup_seqs
  }
  
  # add percent identities for each ortholog - reference ortholog pair
  orthologs <- orthologs %>%
    rowwise() %>%
    mutate(length_ratio = (ortholog_length / ref_ortholog_length) * 100,
           percent_identity = calculate_percent_identity(
             orthogroup_seq_hash[[orthogroup]][[ortholog]],
             orthogroup_seq_hash[[orthogroup]][[ref_ortholog]],
             species_is_outgroup
           ))
}

################################################################################

## Helper functions

get_clades_hash <- function(clades_df) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clades_hash <- new.env()
  for (i in 1 : nrow(clades_df)) {
    clades_hash[[clades_df$species[i]]] = clades_df$clade_broad[i]
  }
  
  return(clades_hash)
}

calculate_percent_identity <- function(seq, ref_seq, species_is_outgroup) {
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  if (!species_is_outgroup) {
    MATRIX = 'BLOSUM45'
  } else {
    MATRIX = 'BLOSUM62'
  }
  
  NW_alignment <- pairwiseAlignment(seq, ref_seq,
                                    type = ALIGNMENT_TYPE,
                                    substitutionMatrix = MATRIX,
                                    gapOpening = GAP_OPEN,
                                    gapExtension = GAP_EXTEND)
  
  return(pid(NW_alignment, type = 'PID2'))
}


plot_lengths_vs_percent_id <- function(orthologs, name, out) {
  # ---------------------------------------------------------------------------
  # Saves a plot of relative ortholog lengths of the species to its yeast
  # orthologs, versus the percent identities between the species-yeast
  # ortholog pairs.
  #
  # Arguments:
  #   orthogroups: dataframe containing ortholog pair lengths and percent ids
  #
  #   name: name of the group of species contained in orthogroups dataframe
  #
  #   out: directory to save the resulting plot in
  #
  # ---------------------------------------------------------------------------
  orthologs <- filter(orthologs, !species_is_outgroup)
  title <- str_c('n = ', nrow(orthologs), ' microsporidia-yeast ortholog pairs\n',
                 sum(orthologs$essential_in_ref), ' yeast-essential pairs, ',
                 nrow(orthologs) - sum(orthologs$essential_in_ref), ' non-essential pairs')
  
  orthologs <- orthologs %>%
    mutate(essential_in_ref = ifelse(essential_in_ref,
                                     'Essential',
                                     'Non-essential'))
  
  ggplot(orthologs, aes(y = percent_identity, x = length_ratio,
                        colour = factor(essential_in_ref))) +
    geom_point(alpha = 0.50, ) +
    geom_rug(alpha = 0.50, aes(color = factor(essential_in_ref))) +
    labs(x = 'Length of species ortholog / length of yeast ortholog',
         y = '% Identity (Identical positions / aligned positions',
         title = title,
         color = 'Ortholog essentiality in yeast') +
    theme_bw() +
    theme(axis.title = element_text(color = 'black', size = 14),
          axis.text = element_text(size = 14),
          title = element_text(color = 'black', size = 14),
          legend.text = element_text(size = 14),
          legend.position = 'bottom',
          legend.justification = 'center')
  
  ggsave(plot = plot, filename = str_c(out, '/', name, '.svg'),
         dpi = 600, units = 'in', width = 9, height = 7.33)
}