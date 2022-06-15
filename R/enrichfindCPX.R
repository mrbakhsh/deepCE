  #' enrichfindCPX
  #' @title Functional Enrichment Analysis for Predicted Modules
  #' @param predcpx A data.frame containingpredicted modules resulted from
  #' \code{\link[deepCE]{get_clusters}},
  #' \code{\link[deepCE]{predPPI_MLP}} and \code{\link[deepCE]{predPPI_1D.CNN}}.
  #' @param threshold Custom p-value threshold for significance.
  #' @param sources A vector of data sources to use.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @param p.corrction.method The algorithm used for multiple testing
  #' correction;defaults to 'bonferroni'.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @param org An organism name;defaults to 'hsapiens'.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @return A data.frame with the enrichment analysis results.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom dplyr arrange
  #' @description This function uses \code{\link[gprofiler2]{gost}} function
  #' in \code{gprofiler2} package to perform functional enrichment analysis
  #' for predicted modules.
  #' @export
  #' @examples
  #' # predict complexes
  #' predcpx <- get_clusters()
  #' predcpx <- predcpx[1:10,]
  #' # perform enrichment for KEGG
  #' enrichCPX <-
  #' enrichfindCPX(predcpx,
  #' sources = c("KEGG"))




  enrichfindCPX <-
    function(predcpx,
             threshold = 0.05,
             sources = c("GO", "KEGG", "CORUM", "REAC", "CORUM"),
             p.corrction.method = "bonferroni",
             org = "hsapiens") {

      colInput <-
        c("ClustID", "Members")
      if(!all(colInput %in% colnames(predcpx))){
        missingCol <-
          setdiff(colInput,
                  colnames(predcpx)[match(colInput, colnames(predcpx))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
      }

      if(!is.data.frame(predcpx)){
        stop("Input data should be data.frame")
      }

      indcpx <-
        str_split(predcpx$Members, "\\s+")
      names(indcpx) <- predcpx$ClustID

      annotCov <- lapply(indcpx, function(x) {
        gprofiler2::gost(x,
                         significant = TRUE,
                         exclude_iea = TRUE,
                         evcodes = TRUE,
                         user_threshold = threshold,
                         sources = sources,
                         correction_method = p.corrction.method,
                         organism = org
        )
      })

      df <-
        lapply(annotCov, function(x) as.data.frame(x[[1]]))
      ans <-
        map_df(df, ~ as.data.frame(.x), .id = "id")

      return(ans)
    }
