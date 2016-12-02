gdina_control <- function(maxit = 2000, crit.min.cng = 1e-6, crit.type = "pj",
                          uniqueMj = TRUE, eps.I = 1e-6, eps.P = 1e-6, trace = 0, lambda = NULL, 
                          useoptim = FALSE, init.pj = NULL, init.pa = NULL, init.random = FALSE, 
                          method = "WLS")
{
  rval <- list(maxit = maxit, crit.min.cng = crit.min.cng, 
               crit.type = crit.type, uniqueMj = uniqueMj, eps.I = eps.I, 
               eps.P = eps.P, trace = trace, lambda = lambda, 
               useoptim = useoptim, init.pj = init.pj, init.pa = init.pa, 
               init.random = init.random, method = method)
  return(rval)
}

gdina <- function(x, q, rule = "G-DINA", link = "identity", 
                  ctrl = gdina_control())
{
  
  link <- match.arg(link, c("identity", "log", "logit"))
  
  linkfn <- switch(link, 
                   "identity" = function(x) x, 
                   "log" = log, 
                   "logit" = qlogis)
  
  invlinkfn <- switch(link, 
                      "identity" = function(x) x, 
                      "log" = exp, 
                      "logit" = plogis)
  
  x <- as.matrix(x) # data
  q <- as.matrix(q) # q-matrix
  
  i <- nrow(x) # number of persons
  j <- nrow(q) # number of items
  k <- ncol(q) # number of skills
  m <- 2^k     # number of attribute patterns in total
  
  prep <- gdina_prepare(q, rule = rule, uniqueMj = ctrl$uniqueMj)
  rule <- prep$rule
  Mj <- prep$Mj
  Kj <- prep$Kj
  a <- prep$a
  np <- prep$np
  alphaeta <- prep$alphaeta
  alphalg <- prep$alphalg
  
  ## Generate unique response patterns and compute their frequency.
  rownames(x) <- apply(x, 1, paste, collapse = "")
  x_patt <- unique(x)
  x_freq <- as.vector(table(rownames(x))[rownames(x_patt)])
  names(x_freq) <- rownames(x_patt)
  
  ## function to compute the conditional likelihood.
  eq14 <- function(pj, xx)
  {
    pja <- t(sapply(1:j, function(jj) pj[[jj]][alphalg[,jj]]))
    lxa <- .Call("calc_lxa", xx, pja, PACKAGE = "Rcdm")
    rownames(lxa) <- rownames(xx)
    colnames(lxa) <- rownames(a)
    return(lxa)
  }
  
  ## function to compute the posterior probability
  posterior <- function(par, individuals = FALSE)
  {
    lxa <- eq14(par$pj, x_patt)
    if(individuals) lxa <- lxa[rownames(x),]
    # lxa <- t(exp(t(log(lxa)) + par$lprior))
    # lxa / rowSums(lxa)
    ll <- colSums(t(lxa) * par$pa)
    t(t(lxa/ll) * par$pa)
  }
  
  ## function to compute the log-likelihood value
  loglik <- function(par)
  {
    # pj <- lapply(1:j, function(jj) {
    #   invlinkfn(prep$Mj[[jj]] %*% par$dj[[jj]])
    # })
    lxa <- eq14(par$pj, x_patt)
    sum(x_freq * log(colSums(t(lxa) * par$pa)))
  }
  
  ## Expectation  
  E_step <- function(par)
  {
    post <- posterior(par) * x_freq
    R <- I <- lapply(2^Kj, function(x) rep(NA, x))
    for(jj in 1:j) {
      rval <- post %*% alphaeta[[jj]]
      I[[jj]] <- colSums(rval)
      R[[jj]] <- colSums(rval * x_patt[,jj])
    }
    return(list(I = I, R = R, post = post))
  }
  
  ## function that performs the E and the M steps
  EM <- function(par, eps = ctrl$eps.P)
  {
    ## expectation
    tval <- E_step(par)
    ## maximization
    pj_upd <- lapply(1:j, function(jj) {
      est <- tval$R[[jj]] / tval$I[[jj]]
      est[is.nan(est)] <- eps
      est
    })
    ## avoid estimates on the boundary
    pj_upd <- lapply(pj_upd, function(pj) {
      pj[pj < eps] <- eps
      pj[pj > 1 - eps] <- 1 - eps
      pj
    })
    ## estimate skill probabilities
    pa_upd <- colSums(tval$post)/i
    ## avoid estimates on the boundary
    pa_upd[pa_upd < eps] <- eps
    pa_upd[pa_upd > 1 - eps] <- 1 - eps
    pa_upd <- pa_upd / sum(pa_upd)
    return(list(pj = pj_upd, pa = pa_upd))
  }
  
  ## function to generate starting values
  startvalues <- function(init = c(0.1, 0.9), random = FALSE, pj_init = NULL, 
                          pa_init = NULL)
  {
    ## item parameters
    if(is.null(pj_init)) {
      if(random) {
        pj_init <- lapply(1:j, function(jj) {
          lgjj <- nrow(prep$Mj[[jj]])
          rval <- cumsum(runif(lgjj, ctrl$eps.P, 1 - ctrl$eps.P))
          if(any(rval > 1)) rval/sum(rval) else rval
        })
      } else {
        pj_init <- lapply(1:j, function(jj) {
          lgjj <- nrow(prep$Mj[[jj]])
          seq(init[1], init[2], length.out = lgjj)
        })
      }
    }
    ## skill parameters
    if(is.null(pa_init)) pa_init <- rep(1/m, m)
    ## names
    names(pj_init) <- rownames(q)
    names(pa_init) <- rownames(a)
    # list(pj = pj_init, pa = pa_init, lprior = log(pa_init))
    list(pj = pj_init, pa = pa_init)
  }
  
  ## compute response probabilities from item pramaters
  pj2dj <- function(pj, I = NULL, method = ctrl$method, lambda = ctrl$lambda)
  {
    rval <- lapply(1:j, function(jj) {
      Mjj <- Mj[[jj]]
      Pjj <- pj[[jj]]
      Wjj <- if(method == "WLS") I[[jj]] else rep(1, length(Pjj))
      if(is.null(lambda)) {
        Gjj <- Hjj <- NULL
        np <- ncol(Mjj)
        if(link == "identity") {
          ## linear predictor is between 0 and 1
          Gjj <- Mjj
          Hjj <- rep(2*ctrl$eps.P, nrow(Mjj))
          Gjj <- rbind(Gjj, -Gjj)
          Hjj <- c(Hjj, Hjj-1)
          ## constraint: intercept between 0 and 1
          Gjj <- rbind(Gjj, c(1L, rep(0L, np-1)))
          Hjj <- c(Hjj, 2*ctrl$eps.P)
          Gjj <- rbind(Gjj, c(-1L, rep(0L, np-1)))
          Hjj <- c(Hjj, 2*ctrl$eps.P-1)
          ## constraint: other parameters between -1 and 1
          Gjj <- rbind(Gjj, diag(rep(1L, np))[-1,])
          Hjj <- c(Hjj, rep(2*ctrl$eps.P-1, np-1))
          Gjj <- rbind(Gjj, diag(rep(-1L, np))[-1,])
          Hjj <- c(Hjj, rep(2*ctrl$eps.P-1, np-1))
        }
        if(link == "log") {
          ## linear predictor is below 0 -> response prob in (0,1)
          Gjj <- -Mjj
          Hjj <- rep(-linkfn(1 - 2*ctrl$eps.P), nrow(Mjj))
          ## constraint: parameters below 1
          Gjj <- rbind(Gjj, diag(rep(-1L, np))[-1,])
          Hjj <- c(Hjj, rep(2*ctrl$eps.P-1, np-1))
        }
        if(is.null(Gjj)) {
          dj_est <- stats::lm.wfit(Mjj, linkfn(Pjj), Wjj, offset = NULL)$coefficients
        } else {
          Ajj <- diag(sqrt(Wjj)) %*% Mjj
          Bjj <- diag(sqrt(Wjj)) %*% linkfn(Pjj)
          dj_est <- limSolve::lsei(A = Ajj, B = Bjj, G = Gjj, H = Hjj, type = 2)$X
        }
        ## lm.fit produces NA values, when the system is (nearly) singular. This
        ## happens, when two estimates are identical. Then, try to compute 
        ## estimates manually:
        naTRUE <- any(is.na(dj_est))
        useoptim <- ctrl$useoptim
        if(naTRUE) useoptim <- TRUE
        if(rule[jj] == "ACDM") useoptim <- TRUE
        if(useoptim) {
          if(naTRUE) {
            dj_init <- startvalues()$dj[[jj]]
          } else {
            dj_init <- dj_est
          }
          Hjj <- Hjj - ctrl$eps.P
          dj_est <- optim_par(Mjj = Mjj, Pjj = Pjj, Wjj = Wjj, ui = Gjj, ci = Hjj, dj_init = dj_init)
        }
      } else {
        out <- glmnet::glmnet(Mjj, Pjj, family = "gaussian", weights = Wjj, lambda = lambda)
        dj_est <- coef(out)[-2,]
      }
      dj_est
    })
    names(rval) <- names(pj)
    rval
  }
  
  optim_par <- function(Mjj, Pjj, Wjj, ui, ci, dj_init)
  {
    Wjj12 <- diag(sqrt(Wjj))
    fn <- function(dj) sum((Wjj12 %*% (linkfn(Pjj) - Mjj %*% dj))^2)
    gradfn <- function(dj) sum(t(Wjj12 %*% (linkfn(Pjj) - Mjj %*% dj)) %*% (Wjj12 %*% Mjj))
    if(is.null(ui)) {
      out <- optim(dj_init, fn, method = "L-BFGS-B", control = list(maxit = 100))
    } else {
      out <- constrOptim(theta = dj_init, f = fn, grad = gradfn, method = "BFGS",
                         ui = ui, ci = ci, control = list(maxit = 100))
    }
    out$par
  }
  
  ## generate staring values
  par_upd <- startvalues(random = ctrl$init.random, pj_init = ctrl$init.pj, 
                         pa_init = ctrl$init.pa)
  
  ## prepare convergence criteria
  diff <- 1
  
  ## prepare for parameter history
  hist <- list(pj = unlist(par_upd$pj), pa = unlist(par_upd$pa))
  
  if(ctrl$trace > 0) cat("Iter\tCriterion\n")
  iter <- 1
  while (iter < ctrl$maxit & diff > ctrl$crit.min.cng)
  {
    
    if(ctrl$trace > 0) cat(iter, "\t", crit, "\n")
    
    ## update
    par <- par_upd
    
    ## EM iteration
    par_upd <- EM(par)
    
    diff <- switch(ctrl$crit.type,
                   "pj" = max(abs(unlist(par$pj) - unlist(par_upd$pj))),
                   "ll" = abs(loglik(par) - loglik(par_upd)) / loglik(par))
    
    ## History
    hist$pj <- rbind(hist$pj, unlist(par_upd$pj))
    hist$pa <- rbind(hist$pa, unlist(par_upd$pa))
    
    ## progress
    iter <- iter + 1
    
  }
  
  if(iter == ctrl$maxit) warning("Maximum number of iterations reached!")
  
  I <- E_step(par)$I
  R <- E_step(par)$R
  
  par_upd$dj <- pj2dj(pj = par_upd$pj, I = I)
  
  rval <- lapply(1:j, function(jj) {
    data.frame(item   = rep(rownames(q)[jj], , np[jj]),
               itemno = rep(jj, np[jj]), 
               rule   = rep(rule[jj], np[jj]), 
               est    = par_upd$dj[[jj]])
  })
  cftable <- do.call("rbind", rval)
  
  ## return object
  rval <- list(call = match.call(),
               x = x, q = q,
               cftable = cftable,
               dj = par_upd$dj,
               pj = par_upd$pj,
               pa = par_upd$pa,
               I = I,
               R = R,
               posterior = posterior(par, individuals = TRUE),
               niter = iter,
               prep = prep,
               npar = sum(np) + m - 1,
               nobs = NROW(x),
               loglik = loglik(par_upd),
               method = ctrl$method,
               parhist = hist,
               link = link,
               ctrl = ctrl)
  
  names(rval$dj) <- rownames(q)
  names(rval$pj) <- rownames(q)
  
  class(rval) <- "gdina"
  return(rval)
  
}

estfun.gdina <- function(object, prob = FALSE, simplify = FALSE)
{
  
  x <- object$x
  q <- object$q
  
  i <- nrow(x)
  j <- nrow(q)
  l <- 2^ncol(q)
  
  dinvlinkfn <- switch(object$link, 
                       "identity" = function(x) rep(1, length(x)), 
                       "log" = function(x) exp(x), 
                       "logit" = function(x) exp(-x) / (1 + exp(-x))^2)
  
  post <- object$posterior
  
  ## score for estimated response probability
  score_pj <- lapply(1:j, function(jj) {
    PaX <- post %*% object$prep$alphaeta[[jj]]
    Pja <- outer(rep(1L, i), object$pj[[jj]])
    xij <- outer(x[,jj], rep(1L, ncol(Pja)))
    PaX * (xij - Pja) / (Pja * (1 - Pja))
  })
  
  ## score for estimated response probability
  score_dj <- lapply(1:j, function(jj) {
    dPja <- as.vector(dinvlinkfn(object$prep$Mj[[jj]] %*% object$dj[[jj]]))
    score_pj[[jj]] %*% (object$prep$Mj[[jj]] * dPja)
  })
  
  lxa_over_l <- t(t(post)/object$pa)
  lxa_over_l[,object$pa < 1e-10] <- 0
  score_pa <- lxa_over_l[,-l] - lxa_over_l[,l]
  
  if(simplify == "array") {
    if(prob)
      return(cbind(do.call("cbind", score_pj), score_pa))
    else
      return(cbind(do.call("cbind", score_dj), score_pa))
  } else {
    return(list(dj = score_dj, pj = score_pj, pa = score_pa))
  }

}

vcov.gdina <- function(object, type = c("full", "partial", "itemwise"), 
                       tol = 1e-10, prob = FALSE)
{
  type <- match.arg(type)
  tval <- estfun.gdina(object)
  scores <- if(prob) tval$pj else tval$dj
  if(type == "itemwise") {
    rval <- lapply(1:length(scores), function(jj) {
      qr.solve(crossprod(scores[[jj]]), tol = tol)
    })
    vcov <- as.matrix(Matrix::bdiag(rval))
  } else {
    scores <- do.call("cbind", scores)
    ok <- rep(TRUE, ncol(scores))
    if(type == "full") {
      scores <- cbind(scores, tval$pa)
      ok <- c(ok, object$pa[seq(ncol(tval$pa))] > 1e-10)
    }
    vcov <- matrix(NA, length(ok), length(ok))
    vcov[ok,ok] <- qr.solve(crossprod(scores[,ok]), tol = tol)
    vcov <- as.matrix(Matrix::bdiag(vcov))
  }
  ## remove identical cols and rows
  return(vcov)
}

coef.gdina <- function(x, type = "all", prob = FALSE) {
  match.arg(type, c("item", "all", "skill"))
  if(type=="skill") {
    ans <- x$pa
  } else {
    ans <- if(prob) do.call("c", x$pj) else unlist(x$dj)
    if(type=="all") {
      ans <- c(ans, x$pa)
    }
  }
  return(ans)
}

confint.gdina <- function(object, alpha = 0.05, prob = FALSE, digits = 4)
{
  cf <- if(prob) unlist(object$pj) else unlist(object$dj)
  se <- sqrt(diag(vcov(object, prob = prob)))[1:length(cf)]
  ans <- data.frame(lower = cf - qnorm(1 - alpha/2) * se, 
             upper = cf + qnorm(1 - alpha/2) * se)
  colnames(ans) <- c(
    sprintf("%3.1f %%", 100*(alpha/2)),
    sprintf("%3.1f %%", 100*(1-alpha/2))
  )
  round(ans, digits)
}

plot.gdina <- function(object, max.plot = 16, plotCI = TRUE, ...)
{
  pj <- object$pj
  j <- length(pj)
  par(mfrow = n2mfrow(min(max.plot, j)))
  if(plotCI) ci <- confint(object, prob = TRUE)
  ans <- sapply(1:j, function(jj) {
    b <- barplot(pj[[jj]], main = paste0("Item", jj), ylim = c(0,1), ...)
    if(plotCI) {
      cijj <- ci[ci$itemno == jj,]
      sapply(1:length(pj[[jj]]), function(p) {
        lines(rep(b[p], 2), c(cijj$lower[p], cijj$upper[p]), col = 2)
        lines(b[p] + c(-1, 1) * 0.01, rep(cijj$lower[p], 2), col = 2)
        lines(b[p] + c(-1, 1) * 0.01, rep(cijj$upper[p], 2), col = 2)
      })
    }
  })
  par(mfrow = c(1,1))
}

barplot.gdina <- function(object, main = "Latent class distribution", ...) {
  barplot(object$pa, las = 2, main = main, ...)
}

gdina_sim <- function(n, q, rule = c("DINA", "DINO", "ACDM", "G-DINA"), 
                      pi = 0.5, dj0 = c(0.2, 0.6), uniqueMj = TRUE)
{
  
  rule <- match.arg(rule)
  
  j <- nrow(q) # number of items
  k <- ncol(q) # number of skills
  m <- 2^k     # number of attribute patterns in total
  
  prep <- gdina_prepare(q, rule = rule, uniqueMj = uniqueMj)
  Mj <- prep$Mj
  Kj <- prep$Kj
  
  ## if not defined, make equal latent class probabilities
  p1 <- rep(pi, ncol(q))
  pa <- apply(t(prep$a) * p1 + (1-t(prep$a)) * (1-p1), 2, prod)
  names(pa) <- rownames(prep$a)
  
  ## sample skill vectors
  alpha <- prep$a[sample(1:m, n, replace = TRUE, prob = pa),]
  
  ## make list if not
  if(!is.list(dj0)) dj0 <- rep(list(dj0), j)

  ## true response probabilities according to de la Torre (2013)
  dj_true <- lapply(1:j, function(jj) {
    npjj <- prep$np[jj]
    init <- rep(0, npjj)
    init[1] <- dj0[[jj]][1]
    init[2:npjj] <- rep(dj0[[jj]][2] / (npjj-1), (npjj-1))
    init
  })
  names(dj_true) <- paste("Item", 1:j, sep = "")
  
  pj_true <- lapply(1:j, function(jj) {
    as.vector(Mj[[jj]] %*% dj_true[[jj]])
  })
  names(pj_true) <- paste("Item", 1:j, sep = "")
  
  ## assign attribute patterns to latent groups
  Lalpha <- prep$alphalg[rownames(alpha),]
  
  ## generate responses
  resp <- sapply(1:j, function(jj) rbinom(rep(1, n), 1, prob = pj_true[[jj]][Lalpha[,jj]]))
  colnames(resp) <- names(Kj)
  rownames(resp) <- apply(resp, 1, paste, collapse = "")
  
  return(list(resp = resp, alpha = alpha, truepar = list(pj = pj_true, dj = dj_true, pa = pa)))
  
}

gdina_prepare <- function(q, rule = "DINA", uniqueMj = TRUE)
{
  
  j <- nrow(q) # number of items
  k <- ncol(q) # number of skills
  
  ## matrix with all possible attribute profiles
  a <- permutations(k)
  
  ## number of skills required per item
  Kj <- rowSums(q)
  
  ## if rule is scalar, make rule to vector
  if(length(rule)==1) rule <- rep(rule, j)
  names(rule) <- names(Kj)
  
  Aj <- lapply(1:j, function(jj) {
    permutations(Kj[jj], col.names = colnames(q)[q[jj,]==1])
  })
  names(Aj) <- names(Kj)
  
  ## The saturated design matrix is the matrix that is used to disentangle
  ## the parameters in the G-DINA model from the conditional probability of
  ## a correct response given the reduced attribute patterns.
  design.matrix <- function(k, rule = "G-DINA") {
    Aj <- permutations(k)
    switch(rule,
           "G-DINA" = {
             ret <- sapply(1:nrow(Aj), function(ll) {
               ret <- sapply(1:k, function(kk) combn(Aj[ll,], kk, FUN = prod), simplify = FALSE)
               do.call("c", ret)
             }, simplify = FALSE)
             ans <- cbind(1, do.call("rbind", ret))
             colnames(ans) <- apply(Aj, 1, function(x) paste(c("d_", which(x == 1)), collapse = ""))
             colnames(ans)[1] <- "d_0"
           },
           "DINA" = {
             ans <- cbind(1, c(rep(0, 2^k-1), 1))
             colnames(ans) <- paste0("d_", c(0, paste0(1:k, collapse = "")))
           },
           "DINO" = {
             ans <- cbind(1, c(0, rep(1, 2^k-1)))
             colnames(ans) <- paste0("d_", c(0, paste0(1:k, collapse = "")))
           },
           "ACDM" = {
             ans <- cbind(1, Aj)
             colnames(ans) <- paste0("d_", 0:k)
           })
    rownames(ans) <- apply(Aj, 1, paste, collapse = "")
    ans
  }
  
  ## The object Mj contains the design matrices for each item.
  Mj <- lapply(1:j, function(jj) design.matrix(Kj[jj], rule = rule[jj]))
  names(Mj) <- names(Kj)
  
  ## number of parameters (per item)
  np <- sapply(Mj, ncol)
  
  ## assign attribute patterns to latent groups
  row.match <- function(x, match, fun = function(x) x)
    apply(x, 1, function(xx) fun(apply(match, 1, function(z) 
      isTRUE(all.equal(z,xx, check.names = FALSE)))))
  
  alphaeta <- lapply(1:j, function(jj) {
    if(uniqueMj) {
      tval <- row.match(a[,q[jj,] == 1, drop = FALSE], Aj[[jj]], fun = which)
      rval <- row.match(unique(Mj[[jj]]), Mj[[jj]][tval,])
    } else {
      rval <- row.match(Aj[[jj]], a[,q[jj,] == 1, drop = FALSE])
    }
    rownames(rval) <- rownames(a)
    rval
  })
  names(alphaeta) <- rownames(q)
  
  alphalg <- sapply(1:j, function(jj) {
    tval <- row.match(a[,q[jj,] == 1, drop = FALSE], Aj[[jj]], fun = which)
    rval <- if(uniqueMj) {
      row.match(Mj[[jj]][tval,], unique(Mj[[jj]]), fun = which)
    } else tval
  })
  rownames(alphalg) <- rownames(a)
  colnames(alphalg) <- rownames(q)
  
  if(uniqueMj) Mj <- lapply(Mj, function(x) unique(x))
  
  return(list(Aj = Aj, Mj = Mj, Kj = Kj, rule = rule, a = a, np = np, 
              alphaeta = alphaeta, alphalg = alphalg))
  
}

logLik.gdina <- function(x, ...) {
  structure(x$loglik, df = x$npar, class = "logLik")
}

item_level_fit <- function(object, red.model = "DINA", prob = FALSE, 
                       method = "none", ...)
{
  
  if(!all(object$prep$rule == "G-DINA"))
    stop("All items in 'm' must be estimated using the G-DINA rule!")
  
  j <- nrow(object$q)
  
  ## restriction matrices  
  R <- lapply(object$prep$Mj, restriction_matrix, rule = red.model, prob = prob)
  
  ## variance covariance
  Vmat <- vcov(object, prob = prob)
  
  ## parameter indices
  rval <- t(sapply(1:j, function(jj) {
    Rj <- R[[jj]]
    Wj <- if(object$prep$Kj[jj] > 1) {
      parind <- c(0, cumsum(object$prep$np))[jj] + 1:object$prep$np[jj]
      tmat <- if(prob) Rj %*% object$pj[[jj]] else Rj %*% object$dj[[jj]]
      Wj <- as.numeric(t(tmat) %*% solve(Rj %*% Vmat[parind,parind] %*% t(Rj)) %*% tmat)
      c(Wj, nrow(Rj), 1 - pchisq(Wj, nrow(Rj)))
    } else c(NA, NA, NA)
  }))
  
  rval <- cbind(object$prep$Kj, rval)
  
  rval[,4] <- p.adjust(rval[,4], method = method)
  
  rownames(rval) <- rownames(object$q)
  colnames(rval) <- c("Kj", "W value", "df", "Pr(>W)")
  
  class(rval) <- "wald.gdina"
  return(rval)
  
}

print.wald.gdina <- function(object)
{
  printCoefmat(object, has.Pvalue = TRUE)
}

restriction_matrix <- function(Mj, rule = "DINA", prob = FALSE)
{
  np <- ncol(Mj)
  k <- log(nrow(Mj), base = 2)
  if(k > 1) {
    
    ## partially adopted from CDM::contraint_matrix
    
    if (rule == "DINA") {
      R <- matrix(0, nrow = 2^k - 2L, ncol = np)
      if(prob) {
        for (vv in 1:nrow(R)) {
          R[vv, vv] <- 1L
          R[vv, vv + 1L] <- -1L
        }
      } else {
        for (vv in 1:nrow(R)) {
          R[vv, vv + 1L] <- 1L
        }
      }
    }
    
    if (rule == "DINO") {
      R <- matrix(0, nrow = 2^k - 2, ncol = np)
      if(prob) {
        for (vv in 1:nrow(R)) {
          R[vv, vv + 1L] <- 1L
          R[vv, vv + 2L] <- -1L
        }
      } else {
        for (vv in 1:nrow(R)) {
          R[vv, ] <- Mj[vv + 2,] - Mj[vv + 1,]
        }
      }
    }
    
    if (rule == "ACDM") {
      R <- matrix(0, nrow = np - (k + 1), ncol = np)
      if(prob) {
        R <- round(MASS::ginv(Mj), 0)[-c(1:(k+1)),,drop = FALSE]
      } else {
        for (vv in 1:nrow(R)) {
          vv1 <- vv + (k + 1)
          R[vv, vv1] <- 1
        }
      }
    }
    
    return(R)
    
  } else return(NULL)
  
}

permutations <- function(k, col.names = NULL, row.names = NULL)
{
  tval <- as.matrix(expand.grid(replicate(k, c(0L, 1L), simplify = FALSE)))
  #tval <- gtools::permutations(2, k, 0:1, repeats.allowed = TRUE)
  rval <- lapply(0L:k, function(kk) {
    tmp <- tval[rowSums(tval) == kk,,drop=FALSE]
    tmp[order(tmp %*% 10^(k:1L), decreasing = TRUE),,drop = FALSE]
  })
  rval <- do.call("rbind", rval)
  rownames(rval) <- if(is.null(row.names)) apply(rval, 1L, paste, collapse = "") else row.names
  colnames(rval) <- if(is.null(col.names)) paste0("K", 1L:k) else col.names
  return(rval)
}

anova.gdina <- function(object, ..., names = NULL)
{
  objects <- list(object, ...)
  title <- "Analysis of Variance Table\n"
  rval <- data.frame("Npar" = sapply(objects, function(m) m$npar),
                     "logLik" = sapply(objects, logLik),
                     "AIC" = sapply(objects, AIC),
                     "BIC" = sapply(objects, function(m) AIC(m, k = log(nrow(m$x)))))
  ord <- order(rval$Npar)
  rval <- rval[ord,]
  rval$Df <- rval$Npar - c(NA, rval$Npar[-length(objects)])
  rval$Deviance <- (-2) * (c(NA, rval$logLik[-length(objects)]) - rval$logLik)
  dfs <- rval$Df
  vals <- rval$Deviance
  rval <- cbind(rval, "Pr(>Chi)" = pchisq(vals, abs(dfs), lower.tail = FALSE))
  if(is.null(names)) {
    rownames(rval) <- paste("m", 1:length(objects), sep = "")
  } else {
    rownames(rval) <- names[ord]
  }
  structure(rval, heading = c(title), class = c("anova", "data.frame"))
}

stepAIC.gdina <- function(object, k = 2, 
                          check.rules = c("DINA", "DINO", "ACDM", "G-DINA"))
{
  rule_upd <- rule <- object$prep$rule
  rval <- matrix(NA, nrow(object$q), length(check.rules))
  colnames(rval) <- check.rules
  rownames(rval) <- rownames(object$q)
  for(jj in which(rowSums(object$q) > 1)) {
    cat("Probing item", jj, "\n")
    rval[jj,] <- sapply(check.rules, function(rr) {
      if(rr != rule[jj]) {
        ruleA <- rule
        ruleA[jj] <- rr
        mA <- update(object, rule = ruleA)
        AIC(mA, k = k)
      } else AIC(object, k = k)
    })
    rule_upd[jj] <- names(which.min(rval[jj,]))
  }
  cat("\n")
  print(rval)
  mnew <- update(object, rule = rule_upd)
  cat("\n")
  cat("Selected rules per item:\n")
  print(mnew$prep$rule)
  return(mnew)
}

update.gdina <- function(object, x = object$x, q = object$q, 
                         rule = object$prep$rule)
{
  return(gdina(x = x, q = q, rule = rule, ctrl = object$ctrl, link = object$link))
}

print.gdina <- function(x)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of items:", nrow(x$q))
  cat("\nNumber of skills:", ncol(x$q))
  cat("\nNumber of latent classes:", 2^ncol(x$q))
  cat("\nNumber of respondents:", x$nobs)
  cat("\nNumber of parameters:", x$np)
  cat("\n\nlog Lik. = ", logLik(x), ", BIC = ", BIC(x), ", AIC = ", AIC(x), sep = "")
}

summary.gdina <- function(x, prob = FALSE)
{
  
  j <- nrow(x$q)
  if(prob) {
    nc <- sapply(x$pj, length)
    rval <- lapply(1:j, function(jj) {
      data.frame(item    = rep(rownames(x$q)[jj], nc[jj]),
                 itemno  = rep(jj, nc[jj]),
                 pattern = rownames(x$prep$Mj[[jj]]),
                 rule    = rep(x$prep$rule[jj], nc[jj]),
                 prob    = x$pj[[jj]])
    })
    cftable <- do.call("rbind", rval)
    cftable$se <- sqrt(diag(vcov(x, prob = TRUE))[1L:nrow(cftable)])
    rownames(cftable) <- NULL
    colnames(cftable) <- c("Item", "Itemno", "Pattern", "Rule", "Estimate", "Std.Err")
  } else {
    np <- sapply(x$dj, length)
    rval <- lapply(1:j, function(jj) {
      data.frame(item   = rep(rownames(x$q)[jj], np[jj]),
                 itemno = rep(jj, np[jj]), 
                 name   = colnames(x$prep$Mj[[jj]]),
                 rule   = rep(x$prep$rule[jj], np[jj]), 
                 est    = x$dj[[jj]])
    })
    cftable <- do.call("rbind", rval)
    cftable$se <- sqrt(diag(vcov(x))[1L:nrow(cftable)])
    rownames(cftable) <- NULL
    colnames(cftable) <- c("Item", "Itemno", "Name", "Rule", "Estimate", "Std.Err")
  }
  
  ans <- list(call = x$call, 
              J = nrow(x$q), 
              K = ncol(x$q), 
              L = 2^ncol(x$q),
              npar = x$npar,
              nobs = x$nobs,
              cftable = cftable,
              pa = x$pa,
              ll = logLik(x),
              aic = AIC(x),
              bic = BIC(x))
  
  class(ans) <- "summary.gdina"
  ans
}

print.summary.gdina <- function(x, digits = 3)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of items:", x$J)
  cat("\nNumber of skills:", x$K)
  cat("\nNumber of latent classes:", x$L)
  cat("\nNumber of respondents:", x$nobs)
  cat("\nNumber of parameters:", x$npar)
  cat("\n\nItem parameters:\n\n")
  print(x$cftable, digits = digits)
  cat("\nSkill parameters:\n\n")
  print(x$pa, digits = digits)
  cat("\nlog Lik. = ", x$ll, ", BIC = ", x$bic, ", AIC = ", x$aic, sep = "")
}

difwald <- function(objR, objF, parm = seq_along(coef(objR, type = "item")), 
                    item = NULL, ...) {
  
  if(any(!objR$prep$rule == "DINA")) stop("Currently, all item must be DINA to perform this test!")
  
  if(!is.null(item)) {
    if(is.numeric(item))
      parm <- which(objR$cftable$itemno %in% item)
    else
      parm <- which(objR$cftable$item %in% item)
  }

  df <- length(parm)
  
  V0 <- matrix(0, df, df)
  betaR <- objR$cftable$est[parm]
  betaF <- objF$cftable$est[parm]
  b <- matrix(c(betaR, betaF))
  R <- cbind(diag(rep(1, df)), diag(rep(-1, df)))
  
  vcovR <- vcov(objR, ...)[parm,parm]
  vcovF <- vcov(objF, ...)[parm,parm]

  # vcovR <- if(!is.null(objR$vcov)) objR$vcov[parm,parm] else vcov(objR, ...)[parm,parm]
  # vcovF <- if(!is.null(objF$vcov)) objF$vcov[parm,parm] else vcov(objF, ...)[parm,parm]
  
  if(any(is.na(vcovR)) || any(is.na(vcovF))) {
    warning("NA values in variance covariance matrix")
    return(c(NA, NA))
  }
  
  V <- rbind(cbind(vcovR, V0), cbind(V0, vcovF))
  W <- t(R %*% b) %*% qr.solve(R %*% V %*% t(R)) %*% (R %*% b)
  pval <- pchisq(drop(W), df = df, lower.tail = FALSE)
  
  rval <- list(statistic = as.numeric(W), 
               parameter = df, 
               p.value = pval, 
               method = "Wald-test for DIF detection",
               data.name = "objR, objF")

  names(rval$statistic) = "Wj"
  names(rval$parameter) = "df"
  
  class(rval) <- "htest"
  
  return(rval)
  
}

difscore <- function(obj, z, parm = seq_along(coef(obj, type = "item")),
                     item = NULL) {
  
  if(!is.null(item)) {
    if(is.numeric(item))
      parm <- which(obj$cftable$itemno %in% item)
    else
      parm <- which(obj$cftable$item %in% item)
  }

  strucchange::sctest(obj, 
         order.by = z, 
         scores = function(x) estfun.gdina(x, prob = FALSE, simplify = "array"),
         parm = parm,
         functional = ifelse(is.factor(z), "LMuo", "dmax"))

}
