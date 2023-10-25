# custom function to transpose while preserving names
transpose_df <- function(microbiome) {
  t_df <- data.table::transpose(microbiome)
  colnames(t_df) <- rownames(microbiome)
  rownames(t_df) <- colnames(microbiome)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(microbiome = t_df)
}