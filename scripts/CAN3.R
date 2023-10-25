CAN3<-
function (method = NULL, M, y = NULL, fraction = 0.5, tm = 20, 
          K = 100, pk = NULL, Nmax = 5, amin = 0, amax = 1, maxIter = 500, 
          popSize = 200, parallel = TRUE, monitor = FALSE) 
{
  source("scripts/EQO_ga2.R")
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
    if (method == "bls_c") {
      seed <- sample(1:nrow(M), size)
      M.train <- M[seed, ]
      y.train <- y[seed]
      M.test <- M[-seed, ]
      y.test <- y[-seed]
      out <- EQO_bls(M.train, y.train, Nmax, K)
      s <- rowSums(cbind(rep(0, nrow(M.test)), M.test[, 
                                                      which(out$x == 1)]))
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
      out <- EQO_ga2("c", M.train, y.train, pk, Nmax, amin, 
                    amax, maxIter, popSize, parallel, monitor)
      print(paste0("Dimensions of M.test is ", dim(M.test)))
      s <- rowSums(cbind(rep(0, nrow(M.test)), M.test[, 
                                                      which(out$x == 1)]))
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
    }
  })
  xv.pileup <- lapply(xv.all, function(x) {
    pairs <- expand.grid(x$taxa, x$taxa)
    pairs$value <- x$perform
    return(pairs)
  })
  xv.all.dt <- data.table::rbindlist(xv.all)
  nodes <- as.data.frame(xv.all.dt[, sum(perform), by = taxa])
  nodes <- nodes[order(nodes$V1, decreasing = TRUE), ]
  colnames(nodes) <- c("taxa", "value")
  nodes$value <- nodes$value/max(nodes$value)
  rownames(nodes) <- nodes$taxa
  nodes$id <- rownames(nodes)
  nodes$label <- rownames(nodes)
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
  return(list(output = xv.all, nodes = nodes, edges = edges))
}
