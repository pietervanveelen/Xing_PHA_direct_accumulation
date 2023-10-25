# note that only the GA option for continuous data is implemented
# the uniform has not been tested yet too.
# note to self, the bls option needs a Gurobi optimizer (commercial package) to work, likely needed for larger datasets

CAN2 <- function (method = NULL, M, y = NULL, fraction = 0.5, tm = 20, 
          K = 100, pk = NULL, Nmax = Nmax, amin = 0, amax = 1, maxIter = 500, 
          popSize = 200, parallel = TRUE, monitor = FALSE, trait_name = NULL, proj=NULL) {
  
  if(missing(proj)) {
    if(!exists("proj", envir = globalenv(), inherits = F)){proj = "ProjID"}
    else {proj = .GlobalEnv$proj}
  }

  if(missing(trait_name)){
    if(!exists("trait_name", envir = parent.frame(), inherits = F)){trait_name = "trait"}
    else {p_frame = parent.frame()
          proj = p_frame$trait_name}
  }
  
  # integer for Nmax (preferably derived from best_Nmax() function)
  if(missing(Nmax)){
    if(!exists("Nmax", envir = parent.frame(), inherits = F)){
      Nmax = 10
      cat("Nmax not provided. Defaults to 10 taxa")
    }
    else {p.frame = parent.frame()
    Nmax = p.frame$Nmax}
  }
  
  if (missing(y)) {
    method <- "ga_u"
  }
  if (missing(pk)) {
    pk <- rep(0, ncol(M))
  }
  M = as.matrix(M)
  size <- round(nrow(M) * fraction, digits = 0)
  
  xv.all <- lapply(1:tm, function(t) {
    print(paste0(paste0(paste0("Running cross-validation ", 
                               t), " out of "), tm))
  tryCatch(
    {  
    if (method == "bls_c") {
      seed <- sample(1:nrow(M), size)
      M.train <- M[seed, ]
      y.train <- y[seed]
      M.test <- M[-seed, ]
      y.test <- y[-seed]
      out <- EQO_bls(M.train, y.train, Nmax, K)
      s <- rowSums(cbind(M.test[,which(out$x == 1)]))
      y.xv <- cor(s, y.test)
      out.xv <- data.table::data.table(taxa = out$members, 
                                       perform = y.xv)
      return(out.xv)
    }
    if (method == "ga_c") {
      seed <- sample(1:nrow(M), size)
      M.train <- M[seed, ]
      y.train <- y[seed]
      M.test <- M[-seed, ]
      y.test <- y[-seed]
      out <- EQO_ga("c", M.train, y.train, pk, Nmax, amin, 
                    amax, maxIter, popSize, parallel, monitor)
      s <- rowSums(M.test[,which(out$x == 1)])      
      y.xv <- cor(s, y.test)
      out.xv <- data.table::data.table(taxa = out$members, 
                                       perform = y.xv)
      return(out.xv)
    }
    if (method == "ga_u") {
      seed <- sample(1:nrow(M), size)
      M.train <- M[seed, ]
      M.test <- M[-seed, ]
      y.train <- 1
      out <- EQO_ga("u", M.train, y.train, pk, Nmax, amin, 
                    amax, maxIter, popSize, parallel, monitor)
      s <- rowSums(cbind(rep(0, nrow(M.test)), M.test[, 
                                                      which(out$x == 1)]))
      y.xv <- mean(s)/sd(s)
      out.xv <- data.table::data.table(taxa = out$members, 
                                       perform = y.xv)
      return(out.xv)
    }},
    error = function(e) {
      cat("Error in iteration", t, ":", conditionMessage(e), "\n")
      return(NULL)}
    )
  })
  
  # remove list elements with NULL in xv.all due to rowsums error, if any
  xv.all_filtered <- Filter(Negate(is.null), xv.all)
  
  # count of effective cross-validations
  n_cross_val = length(xv.all_filtered)
  
  xv.pileup <- lapply(xv.all_filtered, function(x) {
    pairs <- expand.grid(x$taxa, x$taxa)
    pairs$value <- x$perform
    return(pairs)
  })
 
  # compile raw CV data
  result_raw <- 
    map_dfr(xv.all, bind_rows) %>% 
    group_by(perform) %>% 
    mutate(CV_run = cur_group_id()) %>% 
    select(CV_run, everything())
  
  taxa_count = 
    result_raw %>% 
      ungroup() %>% 
      group_by(taxa) %>% 
      dplyr::summarise(taxon_ensemble_prevalence = n(), # how often is a taxon selected as part of the ensemble (of N taxa defined by bestNmax)
                       median_perform = median(perform), # performance = Pearson correlation coefficient between trait and relative abundance of the ensemble the taxon was part of; median is calculated from all cross-validation iterations that the taxon was present in
                       ensemble_preval = 100*(taxon_ensemble_prevalence/n_cross_val)) # an percent index that denotes taxon importance as the prevalence among cross-validations the taxon was part of ensemble.
  
  result_raw <-
    dplyr::left_join(result_raw, taxa_count, by = "taxa") %>% 
    dplyr::mutate(taxa = factor(taxa),
                  taxa = reorder(taxa, -ensemble_preval)) 
    
  result_raw %>%
  write_csv(., file = glue::glue("output_data/{proj}_{trait_name}_CV_data_raw.csv"))
  if(file.exists(glue::glue("output_data/{proj}_{trait_name}_CV_data_raw.csv"))){
    cat("The raw cross-validation data have been saved in output_data/")}
    else {
      stop("Something went wrong with saving raw CV data")
    }
  
  max_CV_n = max(result_raw$CV_run)
  min_r = min(result_raw$perform)

  # plot the relationship for the top performing taxa (i.e. appearing in the ensemble of bestNmax taxa in >10% of CVs)
  plot_raw_cv_data = 
  result_raw %>% 
    mutate(taxa = factor(taxa),
           taxa = reorder(taxa, median_perform)) %>% # reorder taxa by their median pearson correlation across CVs
    filter(ensemble_preval > 10) %>% # filter taxa that are at observed more than 10% of CVs
    ggplot(aes(x=taxa, y=perform)) +
    geom_hline(yintercept = 0, color = "darkred", alpha = 0.5, linetype = "dashed") +
    stat_summary(aes(color = taxa), fun.data = median_hilow, geom = "errorbar", width = 0.2, alpha = 0.3) +
    geom_point(aes(fill = taxa), shape = 21, 
               alpha = 0.5, 
               color = "grey50", 
               position = position_jitterdodge(5)) +
    geom_violin(aes(fill = taxa), alpha = 0.3, color = NA, scale = "area", show.legend = F) +
    stat_summary(aes(color = taxa), fun = mean, geom = "point", shape = 4, size = 3) +
    stat_summary(aes(fill = taxa, size = ensemble_preval), fun = median, geom = "point", shape = 21) +
    labs(x = "Important taxa\n(>10% of CVs)",
         y ="Pearson's r:\ntrait ~ abundance\n[median & 95perc]", 
         title = glue::glue("Ensemble-taxon membership correlations with<br>{trait_name} across {max_CV_n} cross-validations")) +
    scale_fill_manual(values = c(brewer.pal(6, "Set2"), 
                                 brewer.pal(8, "Dark2"), 
                                 brewer.pal(8, "Paired"), 
                                 colorset 
    )) +
    scale_color_manual(values = c(brewer.pal(6, "Set2"), 
                                  brewer.pal(8, "Dark2"), 
                                  brewer.pal(8, "Paired"), 
                                  colorset 
    )) +
    scale_y_continuous(breaks = seq(round(min_r-0.1, digits = 1), 1, 0.1)) +
    scale_size_continuous(name = "Taxon occurrence index\n(%CV in ensemble)") +
    theme_classic() +
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
          axis.title = element_markdown(),
          axis.text.y = element_markdown(),
          plot.title = element_markdown(size = 12),
          legend.position = "bottom") +
    coord_flip() +
    guides(color = "none",
           fill = "none")
  
  print(plot_raw_cv_data)
  # save plot
  ggsave2(plot = plot_raw_cv_data, filename = glue::glue("figures/{proj}_{trait_name}_CV_taxon_importance.pdf"), width = 8, height = 6)
  if(file.exists(glue::glue("figures/{proj}_{trait_name}_CV_taxon_importance.pdf"))){
    cat("Taxon importance plot has been saved successfully")
  }
  
  # create nodes for CAN 
  xv.all.dt <- data.table::rbindlist(xv.all)
  nodes <- as.data.frame(xv.all.dt[, sum(perform), by = taxa])
  nodes <- nodes[order(nodes$V1, decreasing = TRUE), ]
  colnames(nodes) <- c("taxa", "value")
  nodes$value <- nodes$value/max(nodes$value)
  rownames(nodes) <- nodes$taxa
  nodes$id <- rownames(nodes)
  nodes$label <- rownames(nodes)
  
  # create edges for CAN
  pileup.dt <- data.table::rbindlist(xv.pileup)
  weight.dt <- plyr::ddply(pileup.dt, ~Var1 + Var2, dplyr::summarise, 
                           weight = sum(value))
  weight.filter <- weight.dt
  weight.mat <- reshape2::dcast(weight.filter, Var1 ~ Var2, 
                                value.var = "weight")
  rownames(weight.mat) <- as.character(weight.mat$Var1)
  weight.mat <- as.matrix(weight.mat[, -1])
  weight.mat[upper.tri(weight.mat)] <- NA
  diag(weight.mat) <- NA
  edges <- reshape2::melt(weight.mat)
  edges <- edges[!is.na(edges$value), ]
  colnames(edges) <- c("from", "to", "width")
  
  # create final performance overview tibble
  # output performance and importance of taxa


# summarized performances  
final_performance <-  
    xv.all.dt %>% 
    as_tibble() %>% 
    dplyr::group_by(taxa) %>% 
    dplyr::summarise(perform_total = sum(perform),
                     mean_perform = mean(perform),
                     perform_sd = sd(perform),
                     taxon_count_cv_n = n(),
                     taxon_importance = taxon_count_cv_n/n_cross_val*100,
                     perform_var = perform_sd^2) %>% 
    dplyr::mutate(rel_perform = perform_total/max(perform_total)) %>% 
    dplyr::select(taxa, rel_perform, mean_perform, taxon_importance, perform_var, everything()) %>% 
    dplyr::arrange(desc(rel_perform))
  
  
  return(list(output = final_performance, raw_data = result_raw, nodes = nodes, edges = edges, n_cross_val = n_cross_val))
}
