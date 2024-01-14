# -----------------------------------------------------------------------------
#
# Plot ortholog length distributions for outgroup and non-outgroup species
#
# Jason Jiang - Created: Apr/12/2023
#               Last edited: Aug/11/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))

################################################################################

## Global variabLes

OUTGROUPS <- c('R_allo', 'S_pombe', 'D_disc', 'C_eleg', 'H_sapi', 'D_mela', 'D_reri')
EXCLUDED_SP <- c('P_neur', 'D_roes', 'C_dike', 'N_apis', 'N_bomb', 'H_magn', 'H_tvae', 'D_muel')
LENGTH_LEVELS <- c('<30%', '30% - <60%', '60% - <90%', '≥90%')

################################################################################

main <- function() {
  essential_yeast_genes <- read_lines('../../data/essential_yeast_genes.txt')
  
  clades <- get_clades_hash(read_csv('../../data/species_clades.csv'))
  
  orthogroups <- read_csv('../../results/single_copy_orthogroups.csv') %>%
    mutate(essential = ifelse(ref_ortholog %in% essential_yeast_genes,
                              'Essential',
                              'Non-essential'),
           exclude_species = species %in% EXCLUDED_SP,
           length_ratio = ortholog_length / ref_ortholog_length) %>%
    filter(!exclude_species) %>%
    rowwise() %>%
    mutate(clade = clades[[species]])
  
  # plotting for Fig 2B
  length_percents <- orthogroups %>%
    mutate(length_bin = get_length_bin(length_ratio)) %>%
    group_by(clade, essential, length_bin) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(clade) %>%
    mutate(n_clade = sum(n),
           percent_clade = round(n / n_clade * 100, 1))
  
  plot_ortholog_length_distns(length_percents)
  
  # plotting for Fig S2
  # 1) get length bin
  # 2) calculate number of orthologs in each bin for each clade
  # 3) calculate percent within each clade within each bin
  length_percents_supp <- orthogroups %>%
    mutate(length_bin = as.character(cut(length_ratio * 100,
                            breaks = c(0, 10, 20, 30, 40, 50, 60,
                                       70, 80, 90, 100))),
           length_bin = ifelse(is.na(length_bin),
                               ">100%",
                               length_bin)) %>%
    group_by(clade, length_bin) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(clade) %>%
    mutate(n_clade = sum(n),
           percent_clade = round(n / n_clade * 100, 1))
  
  plot_ortholog_length_distns_supp(length_percents_supp)
}

################################################################################

## Helper functions

get_length_bin <- Vectorize(function(length_ratio) {
  if (length_ratio < 0.30) {
    return('<30%')
  } else if (length_ratio >= 0.30 & length_ratio < 0.60) {
    return('30% - <60%')
  } else if (length_ratio >= 0.60 & length_ratio < 0.90) {
    return('60% - <90%')
  } else {
    return('≥90%')
  }
})


get_clades_hash <- function(clades_df) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clades_hash <- new.env()
  for (i in 1 : nrow(clades_df)) {
    clades_hash[[clades_df$species[i]]] = clades_df$clade_broad[i]
  }
  
  return(clades_hash)
}


plot_ortholog_length_distns <- function(length_percents) {
  # plot stacked bar chart of counts of orthologs in each length bin for
  # Microsporidia and outgroups, stacked by essentiality of orthologs
  ggplot(length_percents, aes(x = factor(length_bin, levels = LENGTH_LEVELS),
                            y = percent_clade, fill = essential)) +
    geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
    geom_text(aes(label = str_c(percent_clade, '%')),
              position = position_dodge(width = 0.9),
              vjust = -0.5) +
    labs(x = '% Length to Saccharomyces cerevisiae ortholog',
         y = '% of orthologs within clade',
         fill = 'Essential') +
    scale_fill_manual('Ortholog essentiality in yeast',
                      values = c('#F8766D', '#619CFF')) + # set fill colors
    facet_wrap(~factor(clade, c("Canonical Microsporidia",
                                "Metchnikovellids",
                                "Early Diverging Microsporidia",
                                "Rozella", "Outgroup")), ncol = 3) + # add facet wrap by 'group'
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, color = 'black'), # set color and font size of axis text
      axis.title = element_text(size = 14, color = 'black'), # set color and font size of axis titles
      legend.text = element_text(size = 14, color = 'black'), # set color and font size of legend text
      strip.text = element_text(size = 14, color = 'black'), # set color and font size of facet wrap titles
      legend.justification = 'center',
      legend.position = 'top',
      legend.title = element_text(size = 14),
      title = element_text(color = 'black', size = 14)
    )
}

plot_ortholog_length_distns_supp <- function(length_percents_supp) {
  ggplot(length_percents_supp, aes(x = length_bin,
                                   y = percent_clade,
                                   fill = factor(clade, c("Canonical Microsporidia",
                                                          "Metchnikovellids",
                                                          "Early Diverging Microsporidia",
                                                          "Rozella",
                                                          "Outgroup")))) +
    geom_bar(position = 'dodge2', stat = 'identity', color = 'black') +
    labs(x = '% Length to Saccharomyces cerevisiae ortholog',
         y = '% of orthologs within clade') +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, color = 'black'), # set color and font size of axis text
      axis.title = element_text(size = 14, color = 'black'), # set color and font size of axis titles
      legend.text = element_text(size = 14, color = 'black'), # set color and font size of legend text
      strip.text = element_text(size = 14, color = 'black'), # set color and font size of facet wrap titles
      legend.justification = 'center',
      legend.position = 'top',
      legend.title = element_blank(),
      title = element_text(color = 'black', size = 14)
    )
}

################################################################################

main()