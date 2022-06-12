
  .wtd_rank = function(mat) {
    ranks = apply(mat, 2, rank, ties = "average")
    # weight the ranks
    # calculate the savage scores
    n = nrow(mat)
    reciprocals = 1 / seq_len(n)
    savage = sapply(seq_len(n), function(i) sum(reciprocals[i:n]))
    # replace each rank with the savage score
    savages = ranks
    savages[] = savage[ranks]
    # calculate pearson correlation
    cor = cor(savages, method = "p")
    return(cor)
  }

  .jaccard = function(mat) {
    present = !is.na(mat) & mat > 0
    intersect = crossprod(present)
    union = nrow(mat) - crossprod(!present)
    J = intersect / union
    return(J)
  }

  # Bayes correlation
  .Bayes_Corr <- function(alpha0, beta0, X){
    nrowsX <- nrow(X)
    k <- ncol(X)
    cs <- colSums(X)
    alphas <- matrix(rep(alpha0,nrowsX), nrow=nrowsX, byrow=TRUE) + X
    betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) +
      matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
    alphasPLUSbetas <- alphas + betas

    # First BIG product term for covariance formula
    Psi <- alphas/alphasPLUSbetas -
      matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k,
             byrow=FALSE)

    # Covariance matrix
    cov_mtrx <- Psi %*% t(Psi) / k

    # Variances (this is a column vector of length = nrowsX)
    var_vec <-
      as.matrix( ( rowSums( (alphas*betas)/
                              ( (alphasPLUSbetas^2)*
                                  (alphasPLUSbetas+1) ) )
                   + rowSums(Psi^2) )/k )

    Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
    diag(Bcorrvals) <- 1
    return(Bcorrvals)
  }

  .Bayes_Corr_Prior3 <- function(X){
    d <- dim(X)
    cs <- colSums(X)
    alpha0 <- (cs+1)/(max(cs)+1)
    beta0 <- rep(1,d[2])
    Bcorrvals <- .Bayes_Corr(alpha0,beta0,X)
    return(Bcorrvals)
  }

  .PCCN <- function(x, rept = 10){

    A <- t(x)
    M <-
      ncol(A)
    N <-
      nrow(A)
    i <- 0
    # message("Compute PCC ...")
    pb <-
      txtProgressBar(min = 1, max = rept, style = 3)
    repeat{
      A.rpoisson <-
        apply(A, c(1, 2), function(x){rpois(1, lambda = x)})
      C.rpoisson <-
        A.rpoisson + 1/M
      B.rpoisson <-
        C.rpoisson/rowSums(C.rpoisson)
      i <- i + 1
      setTxtProgressBar(pb, i)
      B.cor <-
        cor(B.rpoisson, use = "pairwise.complete.obs")
      B.cor[is.na(B.cor)] <-0
      if(i == 1){
        PCC.mat <- B.cor
      }else{
        PCC.mat <- PCC.mat + B.cor
      }
      if(i == rept){
        break
      }
    }
    close(pb)
    PCC.mat.avg <-
      PCC.mat/rept
    return(PCC.mat.avg)

  }

  .wcc_f <- function(x,y)
    ptw::wcc(x, y, trwdth = 1)

  #' similarity_score
  #' @title Calculate Pairwise Protein Profile Similarity using Different
  #'   Metrics
  #' @param mat A Co-elution data matrix with proteins in rows and fractions in
  #'   columns.
  #' @param metric  The measure of association including: c(pearson",
  #'   "spearman", "kendall", "bicor", "cosine", "jaccard", "euclidean",
  #'   "manhattan", "weighted_rank", "mutual_information",
  #'   "zero_count_bayesian", "weighted_cross_correlation", or
  #'   "pearson_plusNoise"); Defaults to pearson.
  #' @param rept Poisson iterations for pearson_plusNoise metric, defaults to
  #'   10.
  #' @references Skinnider,M.A. and Foster,L.J. (2021) Meta-analysis defines
  #'   principles for the design and analysis of co-fractionation mass
  #'   spectrometry experiments. Nat. Methods.
  #' @return a matrix of dimensions (# of proteins) x (# of proteins), scoring
  #'   every possible protein pair, in which higher values reflect more similar
  #'   pairs.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom stats cor
  #' @importFrom stats dist
  #' @importFrom utils setTxtProgressBar
  #' @importFrom stats rpois
  #' @importFrom utils txtProgressBar
  #' @importFrom magrittr set_colnames
  #' @description This function computes pairwise protein similarity in a given
  #'  CF-MS experiment using up to 13 indices.
  #' @export
  #' @examples
  #' M1<-matrix(rnorm(36),nrow=6)
  #' M1 <- abs(M1)
  #' rownames(M1) <- c("A","B","C","D","E","F")
  #' scored_Data <- similarity_score(M1, metric = "euclidean")


  similarity_score = function(mat, metric = "pearson", rept = 10) {


    if (!is.matrix(mat)) {
      mat <- as.matrix(mat)
    }

    if(is.character(mat) == TRUE){
      stop("matrix mus include numerical variables")
    }

    if(is.null(row.names(mat))){
      stop("Please specify the row.names")
    }


    # transpose the matrix
    mat = t(mat)

    if (metric %in% c("pearson", "spearman", "kendall")) {
      cor = cor(mat, method = metric, use = 'pairwise.complete.obs')
    } else if (metric == 'bicor') {
      cor = WGCNA::bicor(mat, use = 'pairwise.complete.obs')
    } else if (metric == 'cosine') {
      cor = lsa::cosine(t(mat))
    } else if (metric == 'jaccard') {
      cor = .jaccard(mat)
    } else if (metric == 'euclidean') {
      m <- t(mat)
      scale_c <- apply(m, 2,
                       function(x) (x)/sum(x,  na.rm = TRUE))
      scale_r <- apply(scale_c, 1,
                       function(x) (x)/sum(x, na.rm = TRUE))
      cor =  1- as.matrix(dist(t(scale_r), method = 'euclidean'))
    } else if (metric == 'manhattan') {
      m <- t(mat)
      scale_c <- apply(m, 2,
                       function(x) (x)/sum(x, na.rm = TRUE))
      scale_r <- apply(scale_c, 1,
                       function(x) (x)/sum(x, na.rm = TRUE))
      cor =  1 - as.matrix(dist(t(scale_r), method = 'manhattan'))
    } else if (metric == 'weighted_rank') {
      cor = .wtd_rank(mat)
    } else if (metric == 'mutual_information') {
      cor = WGCNA::mutualInfoAdjacency(mat)$AdjacencyUniversalVersion1
    } else if (metric == 'zero_count_bayesian') {
      cor = .Bayes_Corr_Prior3(t(mat))
    } else if (metric == 'weighted_cross_correlation'){
      cor <-
        1 - as.matrix(proxy::simil(mat,.wcc_f))
    } else if (metric == 'pearson_plusNoise'){
      cor <- .PCCN(mat, rept = rept)
    }  else {
      stop("don't know what to do for metric: ", metric)
    }

    return(cor)
  }

