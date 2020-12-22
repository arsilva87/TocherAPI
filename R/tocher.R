#' Tocher Clustering API
#' 
#' This function computes clusters based on a distance matrix
#' 
#' @param d distance matrix, a dist "object"
#' @param algorithm either "percentual" (default) or "proportion"
#' @author Anderson R da Silva
#' @import biotools
#' @export
#' @examples \dontrun{
#' tocher(d = dist(rnorm(10)))
#' tocher(algorithm = "sequential")
#' }
tocher <- function(d, algorithm = c("original", "sequential"))
  {
    if (!inherits(d, "dist"))
      stop("'d' must be an object of class 'dist'!")
    algorithm <- match.arg(algorithm)
    d <- as.matrix(d)
    n <- nrow(d)
    
    # object labels
    if (is.null(dimnames(d))) {
      lab <- as.character(1:n)
    } else {
      lab <- colnames(d)
    }
    dimnames(d) <- list(lab, lab)
    
    # aux function to find the two closest objects
    fun.min <- function(mat)
    {
      n <- ncol(mat)
      v1 <- v2 <- NULL
      aux <- data.frame(v1 = rep(colnames(mat), each = n),
                        v2 = rep(colnames(mat), times = n),
                        val = as.vector(mat))
      aux2 <- subset(aux, v1 != v2)
      ind <- which.min(aux2[, "val"])
      mi <- aux2[ind, c("v1", "v2")]
      return(c(as.matrix(mi)))
    }
    
    # initial definitions (cluster 1)
    min1 <- fun.min(d)
    g <- list()
    ig <- 1
    g[[ig]] <- min1
    
    # (original) clustering criterion
    d. <- d
    diag(d.) <- NA
    theta <- max(apply(d., 2, min, na.rm = TRUE))
    criterion <- c()
    
    # clustering
    repeat {
      criterion[ig] <- theta
      newlab <- lab[-charmatch(unlist(g), lab)]
      n <- length(newlab)
      if (n < 1) break()
      m <- matrix(0, n + 1, n + 1)
      colnames(m) <- rownames(m) <- c("G", newlab)
      m[newlab, newlab] <- d[newlab, newlab]
      for(j in 1:n) {
        m["G", newlab[j]] <- m[newlab[j], "G"] <-
          mean(d[g[[ig]], newlab[j]])
      }
      comp <- newlab[which.min(m["G", newlab])]
      if (m["G", comp] <= theta) {
        g[[ig]] <- c(g[[ig]], comp)
        # -------------------------------------------
        # forming a new cluster
      } else {
        ig <- ig + 1
        if (n > 1) {
          # theta according to the algorithm
          theta <- ifelse(algorithm == "original", theta,
                          max(apply(d.[newlab, newlab], 2, min, na.rm = TRUE)) )
          newcomp <- fun.min(d[newlab, newlab])
          if (d[newcomp[1], newcomp[2]] <= theta) {
            g[[ig]] <- newcomp
          } else {
            for(i in 1:n) g[[ig + i - 1]] <- newlab[i]
          }
        } else {
          g[[ig]] <- newlab
        }
      }
    }
    
    # output
    ng <- length(g)
    names(g) <- paste("cluster", 1:ng)
    class <- NULL
    for(k in 1:ng) {
      g[[k]] <- noquote(g[[k]])
      for(i in 1:ncol(d)) {
        if (any(lab[i] == g[[k]])) class[i] <- k
      }
    }
    nopc <- sapply(g, length)
    dc <- biotools::distClust(as.dist(d), nopc, unlist(g))
    out <- list(call = match.call(), 
                algorithm = algorithm,
                clusters = g, 
                class = class, 
                criterion = criterion, 
                distClust = dc, 
                d = as.dist(d))
    class(out) <- "tocher"
    return(out)
  }

# ----------------------------------------------
print.tocher <- 
  function (x, ...) 
  {
    cat("\n          Tocher's Clustering \n\n")
    n <- attr(x$d, "Size")
    cat("Call: ")
    print(x$call)
    cat("\nCluster algorithm:", x$algorithm,
        "\nNumber of objects:", n,
        "\nNumber of clusters:", length(x$clusters), "\n")
    if (n < 25) {
      id <- maxmat(x$distClust)
      cat("Most contrasting clusters: ", id[1], " and ", id[2],
          ", with \n   average intercluster distance: ", 
          x$distClust[id[1], id[2]], "\n\n", sep = "")
      print(x$clusters)
    } else {
      cat("Most contrasting clusters (first 4 objects):\n")
      id <- maxmat(x$distClust)
      nopc <- sapply(x$clusters[id], length)
      rest <- max(nopc) - nopc
      dd <- data.frame(1:max(nopc))
      for(i in 1:length(nopc)) {
        dd[, i] <- c(x$clusters[[id[i]]], rep("", rest[id[i]]))
      }
      colnames(dd) <- paste(names(x$clusters[id]), 
                            " (size ", nopc, ")", sep = "")
      print(head(dd, 4))
      cat("... with average intercluster distance:", 
          x$distClust[id[1], id[2]], "\n")
    }
    invisible(x)
  }
