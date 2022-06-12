


  .comp.2.cc.fdr2 <- function (data1, data2,
                              method = "pearson",
                              p.adjust.methods = "local")
  {

    get.lfdr <- function(r) {
      fdr.out <- fdrtool::fdrtool(r, statistic="pvalue")
      fdr.out
    }

    compcorr <- function(n1, r1, n2, r2){
      # Fisher's Z-transformation
      # ad hoc process
      num1a <- which(r1 >= 0.99)
      num2a <- which(r2 >= 0.99)
      r1[num1a] <- 0.99
      r2[num2a] <- 0.99
      num1b <- which(r1 <= -0.99)
      num2b <- which(r2 <= -0.99)
      r1[num1b] <- -0.99
      r2[num2b] <- -0.99
      # z
      z1 <- atanh(r1)
      z2 <- atanh(r2)
      # difference
      dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
      # p-value
      pv <- 2*( 1 - pnorm(abs(dz)) )
      return(list(diff=dz, pval=pv))
    }

    row.names(data1) <- sort(row.names(data1))
    row.names(data2) <- sort(row.names(data2))


    cc1 <- cor(t(data1), method = method)
    cc2 <- cor(t(data2), method = method)
    ccc1 <- as.vector(cc1[lower.tri(cc1)])
    ccc2 <- as.vector(cc2[lower.tri(cc2)])
    n1 <- ncol(data1)
    n2 <- ncol(data2)
    n <- nrow(data1)
    N <- n * (n - 1)/2
    p1 <- rep(1, N)
    p2 <- rep(1, N)
    pdiff <- rep(1, N)
    diff <- rep(1, N)
    mol.names <- rownames(cc1)
    p1 <- cor2.test(n1, ccc1)
    p2 <- cor2.test(n2, ccc2)
    pdiff <- compcorr(n1, ccc1, n2, ccc2)$pval
    diff <- ccc1 - ccc2
    pdiff[(is.na(pdiff)) == TRUE] <- 1
    if (p.adjust.methods == "local") {
      p1.lfdr <- get.lfdr(p1)$lfdr
      p2.lfdr <- get.lfdr(p2)$lfdr
      pdiff.lfdr <-get.lfdr(pdiff)$lfdr
    } else if (p.adjust.methods == "BH" | p.adjust.methods == "bh") {
      p1.lfdr <- p.adjust(p1, method=p.adjust.methods)
      p2.lfdr <- p.adjust(p2, method=p.adjust.methods)
      pdiff.lfdr <- p.adjust(pdiff, method=p.adjust.methods)
    } else {
      p1.lfdr <- rep("not adjusted", N)
      p2.lfdr <- rep("not adjusted", N)
      pdiff.lfdr <- rep("not adjusted", N)
    }

    ## generates combination
    myindex <- which((lower.tri(cc1))==TRUE, arr.ind=TRUE)
    mol.names1 <- mol.names[myindex[,2]]
    mol.names2 <- mol.names[myindex[,1]]

    ## concat
    res <- cbind(mol.names1, mol.names2, ccc1, p1,
                 ccc2, p2,
                 pdiff, diff, p1.lfdr, p2.lfdr, pdiff.lfdr)

    res1 <- data.frame(res)
    res1[,-c(1,2)] <- sapply(res1[, -c(1,2)], as.numeric)

    head <- c("P1", "P2", "r1", "p1",
              "r2", "p2", "p (difference)", "(r1-r2)",
              "lfdr (in cond. 1)", "lfdr (in cond. 2)",
              "lfdr (difference)")
    colnames(res1) <- head
    return(res1)
  }





  #' diffPPI
  #' @title Calculation of Differential Protein-Protein Interactions
  #' @param m1 a numeric matrix with proteins in rows and
  #' fractions in columns
  #' @param m2 a numeric matrix with proteins in rows and
  #' fractions in columns
  #' @param term A list containing known complexes/Go annotation/pathways
  #' provided by the user.
  #' @param p.adjust.methods c("local", holm", "hochberg", "hommel",
  #' "bonferroni", "BH", "BY", "fdr", "none"). Defaults to local.
  #' @param pvalue An integer specifying the p-value cutoff for differential
  #' interactions. Defaults to 0.05.
  #' @param r1_r2 A number specifying the cutoff for the correlation
  #' difference between two conditions. All interactions above and below this
  #' threshold are excluded. Defaults to |0.5|.
  #' @param minSize Minimal size of a gene set to test.
  #' All pathways below the threshold are excluded. Defaults to 5.
  #' See \code{\link[fgsea]{fgsea}} for more details
  #' @param maxSize Maximal size of a gene set to test. All pathways above
  #' the threshold are excluded.
  #' See \code{\link[fgsea]{fgsea}} for more details. Defaults to 500.
  #' @param tpath A character string indicating the path to the project
  #' directory. If the directory is
  #' missing, it will be stored in the Temp directory.
  #' @return Return following data sets:
  #'  \itemize{
  #'  \item{differential_PPI_unfilt} - Unfiltered differential interactions.
  #'  \item{differential_PPI_filtered} - High-confidence differential
  #'  interactions defined by user-defined threshold.
  #'  \item{fgseaRes} - fgsea enrichment terms for high-confidence differential
  #'  PPIs
  #'  \item{vec_rank} - Rank file
  #'  \item{term_list} - term_list
  #'  }.
  #' @description This function first computes the differential correlations
  #' of two correlation metrics via ‘DiffCorr’ R package. Identified DIPs can
  #' then be mapped to complexes/Go annotation/pathways provided by the users
  #' via gene set enrichment analysis (GSEA),
  #' using the ‘fgsea’ R package to determine the biological meaning of DIPs.
  #' @references Fukushima, A. Gene (2013) 518, 209-214
  #' @author Matineh Rahmatbakhsh
  #' @importFrom stats pnorm
  #' @importFrom DiffCorr cor2.test
  #' @importFrom stats p.adjust
  #' @export
  #' @examples
  #' # profile 1
  #' data("m1")
  #' # profile 2
  #' data("m2")
  #' # biological term (here is known complexes)
  #' data("refcpx")
  #' diff_output <-
  #' diffPPI(m1,
  #' m2,
  #' refcpx,
  #' minSize = 2,
  #' maxSize = 10000,
  #' tpath = tempdir())

  diffPPI <-
    function(m1, m2, term,
             p.adjust.methods = "local",
             pvalue = 0.05,
             r1_r2 = 0.5,
             minSize = 5,
             maxSize = 500,
             tpath = tempdir()){

      . <- NULL
      PPI <- NULL
      P1 <- NULL
      P2 <- NULL




      if (!is.matrix(m1)) {
        m1 <- as.matrix(m1)
      }
      if(is.character(m1) == TRUE){
        stop("matrix must include numerical variables")
      }
      if(is.null(row.names(m1))){
        stop("Please specify the row.names for m1 ")
      }

      if (!is.matrix(m2)) {
        m2 <- as.matrix(m2)
      }
      if(is.character(m2) == TRUE){
        stop("matrix must include numerical variables")
      }
      if(is.null(row.names(m2))){
        stop("Please specify the row.names for m2 ")
      }
      if(!is.list(term)){
        stop("Term Must Be List")
      }

      m1 <-
        filter_ConsecutivePep(m1, 2)

      m2 <-
        filter_ConsecutivePep(m2, 2)

      prot <-
        intersect(row.names(m1), row.names(m2))

      m1 <- subset(m1, row.names(m1) %in% prot)
      m2 <- subset(m2, row.names(m2) %in% prot)


      diff_PPI <- .comp.2.cc.fdr2(m1, m2, method = "pearson",
                     p.adjust.methods = p.adjust.methods)

      diff_PPI_filt <-
        diff_PPI %>%
        filter(abs(.[, 8]) > r1_r2 & .[, 11] < pvalue) %>%
        unite(PPI, c(P1, P2), sep = "~", remove =FALSE)

      # input for fgsea
      term_list <-
        lapply(term, function(cpx){
          m <-
            apply(combn(cpx, 2),2, sort)
          s <-
            paste(m[1, ], m[2, ], sep = "~")
          return(s)
        })

      # input for fgsea
      vec_rank <- diff_PPI_filt[, 9]
      names(vec_rank) <- diff_PPI_filt[, 1]
      vec_rank <- rev(sort(vec_rank))


      fgseaRes <- fgsea::fgsea(pathways = term_list,
                               stats    = vec_rank,
                               minSize  = minSize,
                               maxSize  = maxSize)


      fname <- file.path(tpath, "differential_PPI_unfilt.txt")
      write.table(diff_PPI, file = fname,
                  row.names=FALSE, col.names=FALSE,
                  sep="\t", quote=FALSE)


      fname <- file.path(tpath, "differential_PPI_filtered.txt")
      write.table(diff_PPI_filt, file = fname,
                  row.names=FALSE, col.names=FALSE,
                  sep="\t", quote=FALSE)

     output <- list()
     output[["differential_PPI_unfilt"]] <- diff_PPI
     output[["differential_PPI_filtered"]] <- diff_PPI_filt
     output[["fgseaRes"]] <- fgseaRes
     output[["vec_rank"]] <- vec_rank
     output[["term_list"]] <- term_list

     return(output)

     }

