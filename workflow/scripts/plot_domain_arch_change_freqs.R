# -----------------------------------------------------------------------------
#
# Plot frequencies of domain architectural change events in orthologs across
# clades to reference species (yeast)
#
# Jason Jiang - Created: 2022/09/19
#               Last edited: 2023/08/07
#
# Reinke Lab - Microsporidia orthologs
#
# -----------------------------------------------------------------------------

library(tidyverse)

################################################################################

CLADE_LEVELS = c('Canonical Microsporidia', 'Metchnikovellids',
                 'Early Diverging Microsporidia', 'Rozella', 'Outgroup')

DA_LEVELS <- c('Lost domain(s)', 'Gained domain(s)', 'Swapped domain(s)',
               'Conserved domain architecture')

################################################################################

main <- function() {
  essential_ref_genes <- read_lines('../../data/essential_yeast_genes.txt')
  
  clades <- get_clades_hash(read_csv('../../data/species_clades.csv', show_col_types = FALSE))
  
  domain_arch_changes <- read_csv('../../results/aligned_ortholog_domain_architectures.csv') %>%
    filter(!exclude_species,
           is.na(short_enough_for_domain_loss) | short_enough_for_domain_loss,
           !is.na(aligned_domain_archs)) %>%
    rowwise() %>%
    mutate(clade = clades[[species]]) %>%
    ungroup() %>%
    mutate(ref_ortholog_essential = ifelse(ref_ortholog %in% essential_ref_genes,
                                           'Essential',
                                           'Non-essential'),
           DA_conservation = get_DA_conservation(lost_doms, gained_doms, swapped_doms)) %>%
    separate_rows(DA_conservation, sep = '; ') %>%
    select(species, clade, ortholog, ref_ortholog, ref_ortholog_essential,
           aligned_domain_archs, DA_conservation)
}

################################################################################

## Helper functions

get_clades_hash <- function(species_clades) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clades <- new.env()
  map2_chr(species_clades$species,
           species_clades$clade_broad,
           function(x, y) {clades[[x]] <- y})
  
  return(clades)
}


get_DA_conservation <- Vectorize(function(lost_doms, gained_doms, swapped_doms) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (is.na(lost_doms) & is.na(gained_doms) & is.na(swapped_doms)) {
    return("Conserved domain architecture")
  }
  
  DA_changes <- c()
  
  if (!is.na(lost_doms)) {
    DA_changes <- c(DA_changes, 'Lost domain(s)')
  }
  
  if (!is.na(gained_doms)) {
    DA_changes <- c(DA_changes, 'Gained domain(s)')
  }
  
  if (!is.na(swapped_doms)) {
    DA_changes <- c(DA_changes, 'Swapped domain(s)')
  }
  
  return(str_c(DA_changes, collapse = '; '))
}, vectorize.args = c('lost_doms', 'gained_doms', 'swapped_doms'))


plot_DA_change_rates <- function(domain_arch_changes) {
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  clade_counts <- domain_arch_changes %>%
    group_by(clade) %>%
    summarise(n_clade = n())
  
  domain_arch_change_counts <- domain_arch_changes %>%
    group_by(clade, ref_ortholog_essential, DA_conservation) %>%
    summarise(n = n()) %>%
    left_join(clade_counts, by = "clade") %>%
    mutate(percent = round(100 * (n / n_clade), 1))
  
  ggplot(domain_arch_change_counts, aes(x = factor(DA_conservation,
                                                   levels = DA_LEVELS),
                                        y = percent, fill = ref_ortholog_essential)) +
    geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
    ylim(0, 60) +
    geom_text(aes(label = str_c(percent, '%')),
              position = position_dodge(width = 1),
              vjust = -0.5) +
    labs(y = '% of ortholog domain architectural changes within clade',
         fill = 'Essential') +
    scale_fill_manual('Ortholog essential in yeast',
                      values = c('#F8766D', '#619CFF')) + # set fill colors
    facet_wrap(~factor(clade, c("Canonical Microsporidia",
                                "Metchnikovellids",
                                "Early Diverging Microsporidia",
                                "Rozella", "Outgroup")), ncol = 2) + # add facet wrap by 'group'
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, color = 'black'), # set color and font size of axis text
      axis.text.x = element_text(size = 14, color = 'black', angle = 90),
      axis.title = element_text(size = 14, color = 'black'), # set color and font size of axis titles
      axis.title.x = element_blank(),
      legend.text = element_text(size = 14, color = 'black'), # set color and font size of legend text
      strip.text = element_text(size = 14, color = 'black'), # set color and font size of facet wrap titles
      legend.justification = 'center',
      legend.position = 'top',
      legend.title = element_text(size = 14),
      title = element_text(color = 'black', size = 14)
    )
}

################################################################################

main()