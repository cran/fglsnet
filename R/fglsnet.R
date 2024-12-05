#' A Feasible Generalized Least Squares Estimator for Regression Analysis of Outcomes with Network Dependence
#'
#' \code{fglsnet} estimates a multivariate regression model for analyzing outcomes with network dependence. 
#' One nice feature of the function is that it can distinguish three types of error dependence, 
#' including triadic dependence, mutual dependence, and asymmetric dependence.
#'  
#' @param formula A formula indicating the regression model. 
#'
#' @param M The dependence network. 
#' 
#' @param directed Whether the dependence network is directed or undirected.
#' 
#' @param mcorr Whether request multiple correlation coefficients to distinguish triadic, mutual, and asymmetric error dependence.
#' 
#' @param CSE Whether use clustered standard error for the residual regression. Default cluster is the ego unit.
#' 
#' @param k The number of iterations in the fgls estimation.
#' 
#' @param data The data that are used for the regression.
#' 
#' @details The function estimates a multivariate regression model for analyzing outcomes with network dependence.
#' One nice feature of the function is that it can distinguish three types of error dependence, 
#' including triadic dependence, mutual dependence, and asymmetric dependence.
#' 
#' @return A list containing the coefficient \code{coef}, the testing results on the error correlations \code{rtest}, 
#' the estimated error variance \code{Sigma}, the estimated error correlation matrix \code{Omega}, and the OLS estimates \code{ols}.
#' 
#' @export
#' 
#' @import stats 
#' 
#' @import sna
#' 
#' @import network
#' 
#' @import sandwich 
#' 
#' @importFrom Matrix nearPD
#' 
#' @importFrom matrixcalc is.positive.definite
#' 
#' @importFrom MASS mvrnorm
#' 
#' @importFrom lmtest coeftest coefci
#' 
#' @examples
#' data(dat)
#' 
#' g <- fglsnet(Y~ X-1, M = dat$M, directed = TRUE, mcorr = 1, k = 5, data = dat)
#'
#' g$coef
#' 
#' @references
#' 
#' An, Weihua. 2023. ``A Tale of Twin-Dependence: A New Multivariate Regression Model and an FGLS Estimator for Analyzing Outcomes with Network Dependence." \emph{Sociological Methods and Research} 52(4): 1947-1980. 
#' 
#' Greene, William H. (2008). \emph{Econometric Analysis} (6th edition). New Jersey: Pearson Prentice Hall.\cr
#' 


### GLS Estimator
fglsnet <- function (formula, M = NULL, directed = TRUE, mcorr = TRUE, CSE = FALSE, k = 10, data = data){
  
  ## set up a model.frame
  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula = formula, data = data)
  
  Y <- as.numeric(model.response(mf))
  X <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
  
  ### OLS
  g <- lm(Y ~ X-1)
  gout <- cbind(summary(g)$coef, confint(g, level= 0.95))
  res <- g$residuals
  
  ### Get the covariance matrix
  Tigma <- M
  
  ### Get the edge lists.
  if (directed) { Signet <- network(Tigma)
  } else { Signet <- network(Tigma, directed = FALSE)}
  
  edges <- as.matrix(Signet,matrix.type="edgelist")
  
  ### Get mutual and triadic edges
  ## Get the common neighbors
  mTigma <- mutual(Tigma)
  triad <- t(mTigma) %*% mTigma
  
  w <- nrow(edges)
  me <- rep(0, w)
  te <- rep(0, w)
  h <- 1
  for(i in 1:w) {
    if(Tigma[edges[i,2],edges[i,1]] == 1) {
      me[i] <- 1
      if(triad[edges[i,1],edges[i,2]] >= h) te[i] <- 1
    }
  }
  
  # Get different types of edges
  tedges <- edges[te == 1,]
  medges <- edges[me == 1 & te == 0,]
  aedges <- edges[me == 0,]
  saedges <- cbind(aedges[,2],aedges[,1])
  edgesA <- rbind(aedges, saedges)
  
  N <- nrow(M)
  M1 <- matrix(0, N, N)
  diag(M1) <- 1
  
  tM <- M1
  tM[tedges] <- 1
  
  mM <- M1
  mM[medges] <- 1
  
  aM <- M1
  aM[aedges] <- 1
  
  tnet <- network(tM, directed=FALSE)
  edgesT <- as.matrix(tnet,matrix.type="edgelist")
  
  mnet <- network(mM, directed=FALSE)
  edgesM <- as.matrix(mnet,matrix.type="edgelist")
  
  ### Estimates
  j <- 1
  res <- g$residuals

  while (j <= k ) {
    
    ### Create covariance matrix
    tout <- NULL
    mout <- NULL
    aout <- NULL
    rout <- NULL
    
    trho <- 0
    mrho <- 0
    arho <- 0
    
    if (mcorr == 0) {rhos <- c("Error Corr.") 
    } else {
      rhos <- c("Triadic", "Mutual", "Asymmetric")
    }
    
    if (mcorr == 0) {
      tr1 <- res[edges[,1]]
      tr2 <- res[edges[,2]]
      
      # The Wald test
      treg <- lm (scale(tr1)~ scale(tr2)-1)
      tout <- cbind(summary(treg)$coef,confint(treg, level = 0.95))
      tp1 <- treg$df * summary(treg)$r.squared > qchisq(0.95, df=1)
      # Clustered SE
      if (CSE == 1) {
        tout <- cbind(coeftest(treg, vcov.=vcovCL, cluster = tr1), coefci(treg, vcov.=vcovCL, cluster = tr1))
      }
      trho <- tout[1,1]
      rout <- tout
      
      sTigma <- symmetrize(M)
      Omega <- sTigma * trho
      diag(Omega) <- 1
      
    } else {
      ## Estimate triadic correlation
      if (nrow(tedges) >0) {
        tr1 <- res[edgesT[,1]]
        tr2 <- res[edgesT[,2]]
        treg <- lm (scale(tr1)~ scale(tr2)-1)
        tout <- cbind(summary(treg)$coef,confint(treg, level = 0.95))
        # Clustered SE
        if (CSE == 1) {
          tout <- cbind(coeftest(treg, vcov.=vcovCL, cluster = tr1), coefci(treg, vcov.=vcovCL, cluster = tr1))
        }
        trho <- tout[1,1]
      }
      
      ## Estimate remaining mutual dependence
      if (nrow(medges) >0) {
        mr1 <- res[edgesM[,1]]
        mr2 <- res[edgesM[,2]]
        mreg <- lm (scale(mr1)~ scale(mr2)-1)
        mout <- cbind(summary(mreg)$coef,confint(mreg, level = 0.95))
        # Clustered SE
        if (CSE == 1) {
          mout <- cbind(coeftest(mreg, vcov.=vcovCL, cluster = mr1), coefci(mreg, vcov.=vcovCL, cluster = mr1))
        }
        mrho <- mout[1,1]
      }
      
      ## Estimate asymmetric dependence
      if (nrow(aedges) >0) {
        ar1 <- res[aedges[,1]]
        ar2 <- res[aedges[,2]]
        areg <- lm (scale(ar1)~scale(ar2)-1)
        aout <- cbind(summary(areg)$coef,confint(areg, level = 0.95))
        # Clustered SE
        if (CSE == 1) {
          aout <- cbind(coeftest(areg, vcov.=vcovCL, cluster = ar1), coefci(areg, vcov.=vcovCL, cluster = ar1))
        }
        arho <- aout[1,1]
      }
      
      Omega <- M1
      Omega[tedges] <- trho
      Omega <- symmetrize(Omega)
      Omega[medges] <- mrho
      Omega <- symmetrize(Omega)
      Omega[edgesA] <- arho
      diag(Omega) <- 1
      
      rout <- rbind(tout, mout, aout)
    }
    
    # # Test triadic vs. mutual interference
    if (mcorr == 1) {
      rtest32 <- NULL
      rtest21 <- NULL
      rtest31 <- NULL
      
      if (nrow(tedges) >0 & nrow(medges) >0 ) {
        # z <- (trho-mrho)/ sqrt(tout[1,2]^2 + mout[1,2]^2)
        z <- (atanh(trho) - atanh(mrho))/sqrt(1/(length(treg$residuals)-3) + 1/(length(mreg$residuals)-3))
        rtest32 <- 2*pnorm(-abs(z))
      }
      
      if (nrow(medges) >0 & nrow(aedges) >0) {
        # z <- (mrho-arho)/ sqrt(mout[1,2]^2 + aout[1,2]^2)
        z <- (atanh(mrho) - atanh(arho))/sqrt(1/(length(mreg$residuals)-3) + 1/(length(areg$residuals)-3))
        rtest21 <- 2*pnorm(-abs(z))
      }
      
      if (nrow(tedges) >0 & nrow(aedges) >0) {
        # z <- (trho-arho)/ sqrt(tout[1,2]^2 + aout[1,2]^2)
        z <- (atanh(trho) - atanh(arho))/sqrt(1/(length(treg$residuals)-3) + 1/(length(areg$residuals)-3))
        rtest31 <- 2*pnorm(-abs(z))
      }
      
      rtest <- cbind(rtest32, rtest21, rtest31)
      #colnames(rtest) <- c("Triadic vs. Mutual", "Mutual vs. Asymmetric","Triadic vs. Asymmetric" )[1:length(rtest)]
    }
    
    # Sigi <- solve(Omega)
    if (is.positive.definite(Omega)) { Sigi <- chol2inv(chol(Omega))
    } else { Omega <- nearPD(Omega)$mat; Sigi <- chol2inv(chol(Omega)) }
    # Sigi <- chol2inv(chol(Omega))
    
    ### Estimate the model
    x <- model.matrix(g)
    y <- as.matrix(g$model[,1])
    xtxi <- solve(t(x) %*% Sigi %*% x)
    beta <- xtxi %*% t(x) %*% Sigi %*% y
    
    res <- y - x %*% beta
    R2 <- 1- var(res)/var(y)
    
    print(c("***Current Iteration:", j), quote = F)
    j<- j + 1
  }    
  sig <- sqrt(sum(res^2)/g$df)
  
  SE <- sqrt(diag(xtxi))*sig
  tratio <- as.vector(beta /SE)
  pval <- 2 * pt(abs(tratio), df = g$df, lower.tail=FALSE )
  lci <- beta - qt(.975, g$df) * SE
  uci <- beta + qt(.975, g$df) * SE
  coef1 <- cbind(beta, SE, tratio, pval, lci, uci)
  coef2 <- rbind(coef1,rout[,1:6])
  colnames(coef2) <- c("Estimate", "SE", "t", "pval", "lci", "uci")
  
  coef <- as.matrix(coef2)
  rhoname <- rhos[c(!is.null(tout),!is.null(mout),!is.null(aout))]
  rownames(coef) <- c(rownames(gout), rhoname)
  
  Omega2 <- sig^2 * Omega
  
  if (mcorr == 0) {rtest = NULL}
  
  call <- match.call()
  return( list(call = call, formula = formula, coef = coef, rtest = rtest, Sigma = sig, Omega = Omega, ols = gout) )
}

##### Keep mutual ties
mutual <- function (X, type = "matrix", category11 = NA ) {
  Y <- X  + t(X)
  MX <- X
  MX <- ifelse( Y == 2, 1, 0 )
  if ( !is.na(category11) ) { MX <- ifelse( Y == 2, category11, 0 )}
  diag(MX) <- 0
  return(MX)
}

##### Symmetrize a matrix
symmetrize <-  function (X, type = "matrix" ) {
  Y <- X  + t(X)
  MX <- X
  MX <- ifelse( Y >= 1, 1, 0 )
  return(MX)
}
