# -----------------------------------------------------------------------------
#
# Compare domain versus linker amino acid conservation
#
# Jason Jiang - Created: Apr/13/2023
#               Last edited: Sep/16/2023
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(ggsignif))

################################################################################

main <- function() {
  ref_essential_genes <- read_lines('../../data/essential_yeast_genes.txt')
  
  clades <- get_clades_hash(read_csv('../../data/species_clades.csv', show_col_types = FALSE))
  
  domain_linker_length_ratios <- 
    read_csv('../../results/aligned_ortholog_domain_architectures.csv') %>%
    filter(!exclude_species, !is.na(aligned_domain_archs)) %>%
    rowwise() %>%
    mutate(clade = clades[[species]]) %>%
    ungroup() %>%
    mutate(total_domain_length = get_total_domain_length(domain_lengths),
           ref_total_domain_length = get_total_domain_length(ref_domain_lengths),
           linker_length = ortholog_length - total_domain_length,
           ref_linker_length = ref_ortholog_length - ref_total_domain_length,
           'Domains' = total_domain_length / ref_total_domain_length,
           'Linkers' = linker_length / ref_linker_length,
           'Whole Protein' = ortholog_length, ref_ortholog_length,
           essential_in_ref = ifelse(ref_ortholog %in% ref_essential_genes,
                                     'Essential',
                                     'Non-essential')) %>%
    filter(!is.infinite(`Domains`), !is.infinite(`Linkers`),
           !is.nan(`Domains`), !is.nan(`Linkers`))
  
  # p-value calculations
  domain_vs_linker_essential <- wilcox.test(  # p = 0.0033, bonf corrected
    filter(domain_linker_length_ratios, essential_in_ref == 'Essential')[['Domains']],
    filter(domain_linker_length_ratios, essential_in_ref == 'Essential')[['Linkers']],
    paired = TRUE
  )$p.value
  
  domain_vs_linker_nonessential <- wilcox.test(  # 6.49e-17, bonf corrected
    filter(domain_linker_length_ratios, essential_in_ref == 'Non-essential')[['Domains']],
    filter(domain_linker_length_ratios, essential_in_ref == 'Non-essential')[['Linkers']],
    paired = TRUE
  )$p.value
  
  # reformat dataframe for plotting
  domain_linker_length_ratios_plot <- domain_linker_length_ratios %>%
    pivot_longer(cols = c('Domains', 'Linkers', 'Whole Protein'),
                 names_to = 'type',
                 values_to = 'length_ratio') %>%
    mutate(sqrt_length_ratio = sqrt(length_ratio))
  
  # ggplot(data = domain_linker_length_ratios_plot,
  #        aes(x = type, y = sqrt_length_ratio, fill = essential_in_ref)) +
  #   geom_violin() +
  #   scale_fill_manual('Ortholog essential in yeast?',
  #                     values = c('#F8766D', '#619CFF')) +
  #   labs(title = str_c('n = ', nrow(domain_linker_length_ratios),
  #                      ' microsporidia-yeast single-copy ortholog pairs\n')) +
  #   labs(y = 'Length in yeast ortholog / length in microsporidia ortholog (sqrt)') +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(),
  #         axis.text.x = element_text(size = 14, color = 'black', angle = 45,
  #                                    vjust = 1, hjust = 1),
  #         axis.title.y = element_text(size = 14),
  #         axis.text.y = element_text(size = 14, color = 'black'),
  #         legend.text = element_text(size = 14),
  #         legend.title = element_text(size = 14),
  #         legend.position = 'bottom',
  #         legend.justification = 'center',
  #         title = element_text(color = 'black', size = 14))
  
  # plot median domain/linker length ratios between microsporidia and yeast
  domain_linker_length_ratios_med <- domain_linker_length_ratios %>%
    group_by(clade, essential_in_ref) %>%
    summarise(Domain = median(total_domain_length / ref_total_domain_length),
           Linker = median(linker_length / ref_linker_length),
           'Whole Protein' = median(ortholog_length / ref_ortholog_length),
           n = n()) %>%
    pivot_longer(cols = c('Domain', 'Linker', 'Whole Protein'),
                 names_to = 'type',
                 values_to = 'species_to_yeast_ratio')
  
  ggplot(data = domain_linker_length_ratios_med,
         aes(x = type, y = species_to_yeast_ratio, fill = essential_in_ref)) +
    geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
    scale_fill_manual('Yeast-essential ortholog',
                      values = c('#F8766D', '#619CFF')) +
    labs(title = str_c('n = ', nrow(domain_linker_length_ratios),
                       ' single-copy ortholog pairs to yeast')) +
    labs(y = 'Median of ortholog length / yeast ortholog length within clade') +
    ylim(0, 1.0) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color = 'black', angle = 45,
                                     vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14, color = 'black'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'center',
          title = element_text(color = 'black', size = 14)) +
    facet_wrap(~factor(clade, c("Canonical Microsporidia",
                                "Metchnikovellids",
                                "Early Diverging Microsporidia",
                                "Rozella", "Outgroup")),
               ncol=2) +
    theme(strip.text.x = element_text(size = 14))
  
  percent_id = read_csv("percent_identities.csv") %>%
    group_by(essential_in_ref, clade) %>%
    summarise(median_id = median(percent_identity))
  
  ggplot(data = percent_id,
         aes(x = factor(clade, c("Canonical Microsporidia",
                                 "Metchnikovellids",
                                 "Early Diverging Microsporidia",
                                 "Rozella", "Outgroup")),
             y = median_id,
             fill = essential_in_ref)) +
    geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
    scale_fill_manual('Yeast-essential ortholog',
                      values = c('#F8766D', '#619CFF')) +
    labs(title = str_c('n = ', nrow(percent_id),
                       ' single-copy ortholog pairs to yeast')) +
    labs(y = 'Median % identity to yeast within clade') +
    ylim(0, 40) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color = 'black', angle = 45,
                                     vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14, color = 'black'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = 'top',
          legend.justification = 'center',
          title = element_text(color = 'black', size = 14))
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

get_total_domain_length <- Vectorize(function(domain_lengths) {
  if (is.na(domain_lengths)) {
    return(0)
  }
  
  return(sum(as.integer(str_split(domain_lengths, '; ')[[1]])))
})

################################################################################

main()