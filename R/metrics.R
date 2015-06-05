simple.hamming <- function(trueA, estA){
  A.adj <- trueA
  A.adj[A.adj != 0] <- 1
  
  estA.adj <- estA
  estA.adj[estA.adj != 0] <- 1
  
  sum(abs(A.adj - estA.adj))
}

my.shd <- function(m1, m2){
  m1[m1 != 0] <- 1
  m2[m2 != 0] <- 1
  p <- dim(m1)[2]
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}


metrics <- function(trueA, estA){
  trueA <- as.matrix(trueA)
  estA <- as.matrix(estA)
  estA[abs(estA) > 0] <- 1
  trueA[abs(trueA) > 0] <- 1
  diff <- as.vector(estA) -  as.vector(trueA)
  errs <- as.data.frame(table(diff))
  if(is.element(-1, errs$diff)) FN <- errs[errs$diff == "-1", "Freq"] else FN <- 0
  if(is.element(1, errs$diff)) FP <- errs[errs$diff == "1", "Freq"] else FP <- 0
  P <- sum(trueA)
  N <- nrow(trueA)^2 - sum(trueA)
  TP <- P - FN
  TN <- N - FP
  
  precision <- TP/(TP+FP)
  err <- (FP+FN)/(P+N)
  
  TPR <- TP/P
  FPR <- FP/N
  shd <- my.shd(trueA, estA)
  list(shd = shd, err = err, TPR.Recall = TPR, FPR = FPR, precision = precision)
}


#' Computes some performance metrics for estimate of connectiviy matrix A.
#'
#' @param trueA
#' @param est
#' @param thres
#'
#' @return res
metricsThreshold <- function(trueA, est, thres = seq(0.01, 1, by = 0.01)){
  trueA <- as.matrix(trueA)
  est <- as.matrix(est)
  
  res <- matrix(0, nrow = length(thres), ncol = 6)
  for(t in 1:length(thres)){
    est[abs(est) < thres[t]] <- 0
    res[t,] <- c(thres[t], unlist(metrics(trueA, est)))
  }
  res <- as.data.frame(res)
  res <- res[,-3]
  colnames(res) <- c("Threshold", "SHD", "TPR.Recall", "FPR", "Precision")
  res
}

