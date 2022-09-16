
#Stability test without calculation of the p-value


library(sandwich)
library(strucchange) ## for root.matrix()
library(partykit)    ## for mob_beta_suplm (in order to compute p-values)

mob_grow_fluctests_woP <- function(estfun, z, weights, obj = NULL, cluster = NULL,control,
                                   minsize)
{  
  ## set up return values
  m <- NCOL(z)
  pval <- rep.int(NA_real_, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)
  ## variables to test
  mtest <- if(m <= control$mtry) 1L:m else sort(sample(1L:m, control$mtry))
  
  ## estimating functions (dropping zero weight observations)
  process <- as.matrix(estfun)
  ww0 <- (weights > 0)
  process <- process[ww0, , drop = FALSE]
  z <- z[ww0, , drop = FALSE]
  k <- NCOL(process)
  n <- NROW(process)
  nobs <- if(control$caseweights && any(weights != 1L)) sum(weights) else n
  
  ## scale process
  process <- process/sqrt(nobs) ### scaled process HERE (?) 
  vcov <- control$vcov
  if(is.null(obj)) vcov <- "opg"
  if(vcov != "opg") {
    bread <- vcov(obj) * nobs
  }
  if(vcov != "info") {
    ## correct scaling of estfun for variance estimate:
    ## - caseweights=FALSE: weights are integral part of the estfun -> squared in estimate
    ## - caseweights=TRUE: weights are just a factor in variance estimate -> require division by sqrt(weights)
    meat <- if(is.null(cluster)) {
      crossprod(if(control$caseweights) process/sqrt(weights) else process)
    } else {
      ## nclus <- length(unique(cluster)) ## nclus / (nclus - 1L) * 
      crossprod(as.matrix(apply(if(control$caseweights) process/sqrt(weights) else process, 2L, tapply, as.numeric(cluster), sum)))
    }
  }
  J12 <- root.matrix(switch(vcov,
                            "opg" = chol2inv(chol(meat)),
                            "info" = bread, 
                            "sandwich" = bread %*% meat %*% bread
  ))
  process <- t(J12 %*% t(process)) ## NOTE: loses column names
  
  ## select parameters to test
  if(!is.null(control$parm)) {
    if(is.character(control$parm)) colnames(process) <- colnames(estfun)
    process <- process[, control$parm, drop = FALSE]
  }
  k <- NCOL(process)
  
  ## get critical values for supLM statistic
  from <- if(control$trim > 1) control$trim else ceiling(nobs * control$trim)
  from <- max(from, minsize)
  to <- nobs - from
  lambda <- ((nobs - from) * to)/(from * (nobs - to))
  
  beta <- partykit:::mob_beta_suplm
  
  
  
  logp.supLM <- function(x, k, lambda)
  {
    if(k > 40L) {
      ## use Estrella (2003) asymptotic approximation
      logp_estrella2003 <- function(x, k, lambda)
        -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
      ## FIXME: Estrella only works well for large enough x
      ## hence require x > 1.5 * k for Estrella approximation and
      ## use an ad hoc interpolation for larger p-values
      p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
    } else {
      
      ## use Hansen (1997) approximation
      nb <- ncol(beta) - 1L
      tau <- if(lambda < 1) lambda else 1/(1 + sqrt(lambda))
      beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
      #browser()
      dummy <- beta[,(1L:nb)] %*% x^(0:(nb-1))
      dummy <- dummy * (dummy > 0)
      pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
      if(tau == 0.5) {
        p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
      } else if(tau <= 0.01) {
        p <- pp[25L]
      } else if(tau >= 0.49) {
        p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
        ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
        if(!is.finite(p)) p <- mean(c(pp[1L], pchisq(x, k, lower.tail = FALSE, log.p = TRUE)))
      } else {
        taua <- (0.51 - tau) * 50
        tau1 <- floor(taua)
        p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
        ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
        if(!is.finite(p)) p <- mean(pp[tau1 + 0L:1L])
      }
    }
    return(as.vector(p))
  }
  
  ## compute statistic and p-value for each ordering
  for(i in mtest) {
    zi <- z[,i]
    if(length(unique(zi)) < 2L) next
    if(is.factor(zi)) {
      oi <- order(zi)
      proci <- process[oi, , drop = FALSE]
      ifac[i] <- TRUE
      iord <- is.ordered(zi) & (control$ordinal != "chisq")
      
      ## order partitioning variable
      zi <- zi[oi]
      # re-apply factor() added to drop unused levels
      zi <- factor(zi, levels = unique(zi))
      # compute segment weights
      segweights <- if(control$caseweights) tapply(weights[oi], zi, sum) else table(zi)
      segweights <- as.vector(segweights)/nobs
      
      # compute statistic only if at least two levels are left
      if(length(segweights) < 2L) {
        stat[i] <- 0
        pval[i] <- NA_real_
      } else if(iord) {
        proci <- apply(proci, 2L, cumsum)
        tt0 <- head(cumsum(table(zi)), -1L)
        tt <- head(cumsum(segweights), -1L)
        if(control$ordinal == "max") {  ### ordinal case
          stat[i] <- max(abs(proci[tt0, ] / sqrt(tt * (1-tt))))
          pval[i] <- log(as.numeric(1 - mvtnorm::pmvnorm(
            lower = -stat[i], upper = stat[i],
            mean = rep(0, length(tt)),
            sigma = outer(tt, tt, function(x, y)
              sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
          )^k))
        } else {
          proci <- rowSums(proci^2)
          stat[i] <- max(proci[tt0] / (tt * (1-tt)))
          #pval[i] <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = control$nrep)$computePval(stat[i], nproc = k))
        }
      } else {      
        stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[,j], zi, sum)^2)/segweights))
        #pval[i] <- pchisq(stat[i], k*(length(levels(zi))-1), log.p = TRUE, lower.tail = FALSE)
      }
    } else {
      oi <- if(control$breakties) {
        mm <- sort(unique(zi))
        mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
        order(zi + runif(length(zi), min = -mm, max = +mm))
      } else {
        order(zi)
      }
      proci <- process[oi, , drop = FALSE]
      proci <- apply(proci, 2L, cumsum)
      tt0 <- if(control$caseweights && any(weights != 1L)) cumsum(weights[oi]) else 1:n
      
      from_to <- tt0 >= from & tt0 <= to
      stat[i] <- if(sum(from_to) > 0L) {
        
        xx <- rowSums(proci^2)
        xx <- xx[from_to]
        tt <- tt0[from_to]/nobs
        
        max(xx/(tt * (1 - tt))) 
        
      } else {
        0
      }
      #pval[i] <- if(sum(from_to) > 0L) logp.supLM(stat[i], k, lambda) else NA
    }
  }
  
  ## select variable with minimal p-value
  #best <- which.min(pval)
  
  return(stat)
}


# 
# ## GDP of a logistic regression.
# set.seed(77777L)
# N      <- 10
# b_low  <- 0.1
# b_high <- 25
# x1     <-  rnorm(n = N, mean =  1.5, sd = 4)
# z1     <-  rnorm(n = N, mean =  0, sd = 1)
# z2     <-  rnorm(n = N, mean =  1.5, sd = 4)
# y_lin  <-  x1* ifelse(z1>0 , b_high, b_low)   
# pr     <-  1/(1+exp(-y_lin)) 
# y      <-  rbinom(N,1,pr)  
# clust  <- sample.int(n = 100, size = N, replace = TRUE)
# df <- data.frame(x1 = x1, z1 = z1, z2 = z2, y = y, clust = clust)
# ## Estimated model
# m <- glm( y~x1, data=df, family="binomial")
# ## fluctuation test using sandwich matrix   
# #       control$vcov <- "sandwich"
# 
# ## Inputs for the m-fluctuation test.
# estfun <- sandwich::estfun(m)
# z <- data.frame(z1=z1, z2=z2)
# obj <- m
# weights <- rep(1,nrow(df))
# 
# ## mob$control minimal required options 
# control <- list(vcov =NULL, 
#                 caseweights = F, 
#                 mtry =  Inf, 
#                 trim =  0.1, 
#                 breakties = FALSE )
# ## minsize of leaf
# minsize <-  2
# ## fluctuation test using outer product gradients 
# control$vcov <- "opg"
# mob_grow_fluctests(estfun= estfun, 
#                    z = z, 
#                    weights=weights,
#                    obj = obj, 
#                    cluster = clust)
# 
# # $pval
# #        z1        z2 
# # 0.7224722 0.6902683 
# # 
# # $stat
# #       z1       z2 
# # 3.653720 3.836886 
# # 
# # $best
# # z2 
# # 2 
