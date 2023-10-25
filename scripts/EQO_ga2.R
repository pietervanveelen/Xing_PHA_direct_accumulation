EQO_ga2<-
function (pattern = NULL, M, y = NULL, pk = NULL, Nmax = NULL, 
          amin = 0, amax = 1, maxIter = 500, popSize = 200, parallel = TRUE, 
          monitor = plot) 
{
  if (missing(y)) {
    y <- 1
    pattern <- "u"
  }
  if (missing(pk)) {
    pk <- rep(0, ncol(M))
  }
  if (missing(Nmax)) {
    Nmax <- ncol(M)
  }
  M = as.matrix(M)
  an <- colMeans(M)
  m <- nrow(M)
  n <- ncol(M)
  pen <- sqrt(.Machine$double.xmax)
  c0 <- function(x) {
    ifelse(min(x - pk) < 0, 1, 0)
  }
  c1 <- function(x) {
    (rep(1, n) %*% x) - Nmax
  }
  c2 <- function(x) {
    amin - (an %*% x)
  }
  c3 <- function(x) {
    (an %*% x) - amax
  }
  if (pattern == "c") {
    M0 <- apply(M, 2, function(x) {
      x - mean(x)
    })
    y0 <- y - mean(y)
    P <- t(M0) %*% M0
    Q <- t(M0) %*% y0 %*% t(y0) %*% M0
    Q2 <- t(M0) %*% y0
    fitness <- function(x) {
      (t(x) %*% Q2)/sqrt((t(x) %*% P) %*% x) - (pen * max(c1(x), 
                                                          0)) - (pen * c0(x))
    }
  }
  if (pattern == "d") {
    M0 <- apply(M, 2, function(x) {
      x - mean(x)
    })
    y0 <- y
    L <- matrix(0, nrow = ncol(y0), ncol = ncol(y0))
    diag(L) <- 1/sqrt(colSums(y0))
    P <- t(M0) %*% M0
    Q <- t(M0) %*% y0 %*% L %*% L %*% t(y0) %*% M0
    fitness <- function(x) {
      ((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x) - (pen * 
                                                     max(c1(x), 0)) - (pen * c0(x))
    }
  }
  if (pattern == "u") {
    M0 <- M
    e <- as.matrix(rep(1, m), ncol = 1)
    P <- t(M0) %*% M0 - ((2/n) * (t(M0) %*% e %*% t(e) %*% 
                                    M0)) + ((1/(n^2)) * (t(M0) %*% e %*% t(e) %*% e %*% 
                                                           t(e) %*% M0))
    Q <- ((1/(n^2)) * t(M0) %*% e %*% t(e) %*% M0)
    fitness <- function(x) {
      ((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x) - (pen * 
                                                     max(c1(x), 0)) - (pen * max(c2(x), 0)) - (pen * 
                                                                                                 max(c3(x), 0)) - (pen * c0(x))
    }
  }
  GA <- GA::ga("binary", fitness = fitness, nBits = n, maxiter = maxIter, 
               popSize = popSize, names = colnames(M0), monitor = monitor, 
               parallel = parallel)
  if (pattern == "c") {
    x <- as.numeric(GA@solution)
    fitness <- max(GA@fitness[!is.na(GA@fitness)])
    members <- colnames(M0)[GA@solution == 1]
    abundance <- rowSums(cbind(rep(0, m), rep(0, m), M[, 
                                                       GA@solution == 1]))
    print(paste0("Length of GA@solution is ", ncol(M[,GA@solution == 1])))
    r <- cor(abundance, y0)
    return(list(fitness = fitness, x = x, members = members, 
                abundance = abundance, performance = r))
  }
  if (pattern == "d") {
    solution1 = GA@solution[1, ]
    x <- as.numeric(solution1)
    fitness <- max(GA@fitness[!is.na(GA@fitness)])
    members <- colnames(M)[solution1 == 1]
    abundance <- rowSums(cbind(rep(0, m), rep(0, m), M[, 
                                                       solution1 == 1]))
    s <- abundance - mean(abundance)
    R2 <- (t(s) %*% y0 %*% L %*% L %*% t(y0) %*% s)/(t(s) %*% 
                                                       s)
    return(list(fitness = fitness, x = x, members = members, 
                abundance = abundance, performance = R2))
  }
  if (pattern == "u") {
    x <- as.numeric(GA@solution)
    fitness <- max(GA@fitness[!is.na(GA@fitness)])
    members <- colnames(M0)[GA@solution == 1]
    abundance <- rowSums(cbind(rep(0, m), rep(0, m), M0[, 
                                                        GA@solution == 1]))
    CV <- sd(abundance)/mean(abundance)
    return(list(fitness = fitness, x = x, members = members, 
                abundance = abundance, performance = CV))
  }
}