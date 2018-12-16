#' Learn a Tweedie regression model parameters
#'
#' This is an implementation of the maximum likelihood method
#' to estimate the parameters of Tweedie regression model.
#'
#' @param data a data frame containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.
#' @param pgrid a numeric vector containing the values among which the power
#' parameter p is chosen in order to maximize the model's log-likelihood.
#' If the power parameter is knonw, pgrid will contain a single value.
#'
#' @return This function returns a named list consisting of the Tweedie
#' regression parameters ("p" : the power parameter, "phi": the dispersion
#' parameter, "beta" : the regression coefficients and "glm": a regression
#' object as returned by the glm or lm functions.)
#'
#' @examples
#' \dontrun{
#' learn.trm(data, V5~V2 + V3, pgrid=seq(2,3,0.05))
#' }
#'
#'
#' @export

learn.trm <- function(data, formula, pgrid) {
  max.ll = -Inf
  trm.params = list()
  mf <- model.frame(formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  for (p in pgrid) {
    if (p == 0) {
      res = lm(data = data, formula = formula)
    }
    else
    {
      mod = try({
        glm(
          data = data,
          formula = formula,
          family = tweedie(link.power = 0, var.power = p)
        )
      }, silent = TRUE)
      if (isTRUE(class(mod) == "try-error")) {
        res = NA
        next
      } else{
        res = mod
      }
    }
    mu <- res$fitted.values
    q = ncol(X)
    n = nrow(X)
    phi = sum(((Y - mu) ^ 2) / (mu ^ p)) / (n - q)
    current.ll = scores.trm(
      beta = res$coefficients,
      phi = phi,
      p = p,
      X = X,
      Y = Y
    )[["ll"]]

    if (current.ll > max.ll) {
      max.ll = current.ll
      trm.params[["p"]] = p
      trm.params[["phi"]] = phi
      trm.params[["beta"]] = res$coefficients
      trm.params[["ll"]] = max.ll
      trm.params[["glm"]] = res
    }
  }
  trm.params
}
#' Learn a Tweedie Bayesian Network.
#'
#' This function builds a Tweedie Bayesian network based on a data set.
#' It combines both structure and parameters learning. If a known graph is
#' provided, only the parameters will be estimated relying on the provided
#' network structure.
#'
#' @param data a data frame containing the sample used to learn the model.
#' @param pgrids a list of numeric vectors containing the values among which
#' the power parameters are chosen in order to maximize the model's
#' log-likelihood. If the power parameter is knonw, each pgrid vector will
#' contain a single value.
#' @param G (optional) the Tweedie Bayesian network structure as an adjacency
#' matrix i.e., G[i,j] =1 if Xj is a parent of Xi. If no graph is provided, the
#' structure will be learned relying on student statistical tests.
#' @param nodes.order (optional) The nodes ordering defining an ancestral
#' ordering of the Bayesian network.
#' @param conf.level (optioal) The confidence level for the structure learning.
#' Default value is 0.05.
#'
#' @return This function returns a named list consting of the Tweedie bayesian
#' network structure and parameters.
#'
#' @examples
#' # 3 nodes: the third is node has a Gaussian conditional distribution.
#' \dontrun{
#' learn.tbn(data, pgrids=list(c(2,3), c(2.2,2.3,2.4), c(0)))
#' }
#'
#'
#' @export
learn.tbn <- function(data,
                      pgrids,
                      G = NULL,
                      nodes.order = NULL,
                      conf.level = 0.05) {
  if (!is.null(nodes.order)) {
    data = data[, nodes.order]
  }
  nodes = names(data)
  d = ncol(data)
  n = nrow(data)
  ok = TRUE
  beta = matrix(0, d, d)
  tvalue = beta
  rawbeta = beta
  phi = vector(mode = 'numeric', length = d)
  p   = vector(mode = 'numeric', length = d)
  TBN = list()
  if (is.null(G)) {
    G = matrix(0, d, d)
    #learn structure
    for (i in 2:d) {
      fm = paste(nodes[i], "~", sep = "")
      for (j in 1:(i - 1)) {
        fm = paste(fm, nodes[j], sep = "+")
      }
      trm.params = learn.trm(data = data,
                             formula = as.formula(fm),
                             pgrid = pgrids[[i]])
      res = trm.params$glm
      c = summary(res)
      rawbeta [i, 1:i] = res$coefficients
      G[i, 1:(i - 1)] = as.numeric(c$coefficients[2:i, 4] < conf.level)
    }
  }
  TBN[["Graph"]] = G
  # Estimate parameters using the structure G.
  PaX = matrix(1, n, 1)
  PaX = cbind(PaX, as.matrix(data))
  for (i in 1:d) {
    hasparents = sum(G[i,]) > 0
    if (hasparents) {
      fm = paste(nodes[i], "~", sep = "")
      for (node in nodes[G[i,] > 0]) {
        fm = paste(fm, node, sep = "+")
      }
    }
    else
    {
      fm = paste(nodes[i], "~ 1", sep = "")
    }
    trm.params = learn.trm(data = data,
                           formula = as.formula(fm),
                           pgrid = pgrids[[i]])
    res = trm.params[["glm"]]
    beta[i, c(1, (2:i)[G[i,] > 0])] <- res$coefficients
    tvalue[i, c(1, (2:i)[G[i,] > 0])] <-
      summary(res)$coefficients[, 3]
    phi[i] = trm.params$phi
    p[i] = trm.params$p
  }

  TBN[["beta"]] = beta
  TBN[["rawbeta"]] = rawbeta
  TBN[["phi"]] = phi
  TBN[["power"]] = p
  TBN[["converged"]] = ok
  TBN[["tvalue"]] = tvalue
  TBN
}
#' Compute likelihood-based score for a Tweedie regression model.
#'
#' This function computes the log-likelihood of a Tweedie regression model given
#' a data set. It also provides the MSE, RMSE and NRMSE.
#'
#' @param beta a vector of regression coefficients.
#' @param phi the dispersion parameter of the Tweedie regression model.
#' @param p the power parameter of the Tweedie regression model.
#' @param Y a vector of responses.
#' @param X a matrix of covariates (including a vector of ones).
#'
#' @return This function returns a named list consisting of the Tweedie
#' regression model scores: Mean Squared Error (MSE), Root Mean Squared Error
#' (RMSE), Normalized Root Mean Squared Error (NRMSE) and the log-likelihood
#' (ll).
#'
#' @export
scores.trm <- function(beta, phi, p, Y, X) {
  if (p == 0) {
    mu <- as.numeric(X %*% beta)
    d <- dnorm(x = Y,
               mean = mu,
               sd = sqrt(phi))
  }
  else{
    mu <- as.numeric(exp(X %*% beta))
    d = dtweedie(
      y = Y,
      mu = mu,
      phi = phi,
      power = p
    )
  }
  ll <- sum(log(d[d > 0]))
  mse = mean((Y[d > 0] - mu[d > 0]) ^ 2)
  list(
    MSE = mse,
    RMSE = sqrt(mse),
    NRMSE = sqrt(mse) / (mean(mu[d > 0])),
    ll = ll
  )
}

#' Compute the log-likelihood within a Tweedie Bayesian network.
#'
#' This function computes the log-likelihood within a Tweedie Bayesian network
#' for a given data set.
#'
#' @param beta a matrix of the Tweedie Bayesian network coefficients.
#' @param phi a vector of dispersion parameters.
#' @param p a vector of power parameters.
#' @param X a data frame containing the variables in the model.
#'
#' @return This function returns a named list consisting of the Tweedie
#' Bayesian network scores: Mean Squared Error (MSE), Root Mean Squared Error
#' (RMSE), Normalized Root Mean Squared Error (NRMSE) and the log-likelihood
#' (ll).
#'
#' @examples
#' \dontrun{
#' ll.tbn(beta, phi, p, X)
#' }
#'
#'
#' @export

ll.tbn <- function(beta, phi, p, X) {
  d <- ncol(X) - 1
  res = list()
  MSE = c()
  RMSE = c()
  NRMSE = c()
  logl = c()
  ll = 0
  for (i in 1:d) {
    scores = scores.trm(beta[i, 1:i], phi[i], p[i], X[, i + 1], as.matrix(X[, 1:i]))
    ll = ll + scores[["ll"]]
    MSE = c(MSE, scores[["MSE"]])
    RMSE = c(RMSE, scores[["RMSE"]])
    NRMSE = c(NRMSE, scores[["NRMSE"]])
    logl = c(logl, scores[["ll"]])
  }
  res[["ll"]] = ll
  res[["locScore"]] = data.frame(
    MSE = MSE,
    RMSE = RMSE,
    NRMSE = NRMSE,
    ll = logl
  )
  res
}
#' Computes several likelihood-based score for a given Tweedie Bayesian Network.
#'
#' Computes several likelihood-based score for a given Tweedie Bayesian Network:
#' the log-likelihood, the Akaike Information Criterion (AIC) and the Bayesian
#' Information Criterion (BIC).
#'
#'
#' @param beta a matrix of the Tweedie Bayesian network coefficients.
#' @param phi a vector of dispersion parameters.
#' @param p a vector of power parameters.
#' @param X a data frame containing the variables in the model.
#'
#' @return This function returns a named list consisting of the Tweedie
#' Bayesian network scores: Akaike Information Criterion (AIC), Bayesian
#' Information Criterion (BIC) and the log-likelihood (ll).
#'
#' @export
scores.tbn <- function(beta, phi, p, X) {
  n = nrow(X)
  d = length(p)
  k = sum(beta > 0) + 2 * d
  locll = ll.tbn(beta, phi, p, X)
  ll = locll[["ll"]]
  bic <- log(n) * k - 2 * ll
  aic <- 2 * k - 2 * ll
  list(AIC = aic,
       BIC = bic,
       locScore = locll$locScore)
}
#' Compute the sensitivity of TRM to coefficients modification.
#'
#' Compute the sensitivity of TRM to coefficients modification. The returned
#' value corresponds to the local sensitivity of a Tweedie Bayesian network.
#'
#' @param beta a vector of regression coefficients.
#' @param phi the dispersion parameter of the Tweedie regression model.
#' @param p the power parameter of the Tweedie regression model.
#' @param X a data frame containing the variables in the model.
#' @param db a vector of modification of the Tweedie regression coefficients.
#'
#' @return This function returns the sensitivity of the tweedie regression model
#' induced by the regression coefficients modification.
#'
#' @examples
#' \dontrun{
#' sensi.trm(c(0.2,0,0,1,-2),1,2, c(-0.1,0,0,-0.2,0), X)
#' }
#'
#'
#' @export
sensi.trm <- function(beta, phi, p, db, X) {
  AverageSensi = 0
  if (sum(abs(db)) != 0) {
    if (p == 0) {
      mu1 <- as.numeric(X %*% beta)
      mu2 <- as.numeric(X %*% (beta + db))
    }
    else{
      mu1 <- as.numeric(exp(X %*% beta))
      mu2 <- as.numeric(exp(X %*% (beta + db)))
    }
    if (p == 2) {
      sensi = mu1 / (phi * (1 - p)) * (mu1 ^ (1 - p) - mu2 ^ (1 - p)) - (log(mu1) - log(mu2)) / phi
    } else{
      sensi = mu1 / (phi * (1 - p)) * (mu1 ^ (1 - p) - mu2 ^ (1 - p)) - 1 / (phi *(2 - p)) * (mu1 ^ (2 - p) - mu2 ^ (2 - p))
    }
    AverageSensi = mean(sensi)
  }

  AverageSensi
}

#' Compute a tbn Sensitivity to coefficients modification
#'
#' Computes the sensitity of a Tweedie Bayesian network induced by a
#' modification of the coefficient matrix. This sensistivity relies on the
#' Kulback-Liebler divergence between the joint probability distributions
#' induced by the original and the modified Tweedie bayesian networks.
#'
#' @param beta a lower triangular matrix of the Tweedie Bayesian network coefficients.
#' @param phi a vector of dispersion parameters.
#' @param p a vector of power parameters.
#' @param db an upper-triangular matrix of coefficients modification.
#' @param X a data frame containing the variables of the model.
#'
#' @return This function returns a double corresponding to the Tweedie Bayesian
#' network sensitivity.
#'
#' @examples
#' \dontrun{
#' beta = matrix(0,ncol=5,nrow=5)
#' beta[1,1] = 1;
#' beta[2,1] = 0.5
#' beta[3,1] = -0.5; beta[3,2] = -0.7;beta[3,3] = -1.2
#' beta[4,1] = -2; beta[4,3] = -1
#' beta[5,1] = 0.2; beta[5,4] = 1; beta[5,5] = -2
#' db = matrix(0,ncol=5,nrow=5)
#' db[2,2] = 0.1 # adding an arc from X1 to X2.
#' sensi.tbn(beta, c(1,2,1.5,1,1), c(2,2,3,2,0), db, X)
#' }
#'
#' @export

sensi.tbn <- function(beta, phi, p, db, X) {
  d <- ncol(X) - 1
  res = list()
  s = 0
  localSensis = vector(mode = "numeric", length = nrow(X))
  for (i in 1:d) {
    temp.res = sensi.trm(beta[i, 1:i],
                         phi[i],
                         p[i],
                         db[i, 1:i],
                         as.matrix(X[, 1:i]))
    s = s + temp.res
  }
  s
}
#' Generates a sample from Tweedie Bayesian network
#'
#' Generates an n-random sample according to a Tweedie Bayesian network with
#' parameters (p,phi,beta)
#'
#' @param beta a lower triangular matrix of the Tweedie Bayesian network coefficients.
#' @param phi a vector of dispersion parameters.
#' @param p a vector of power parameters.
#' @param n an integer specifying the sample size.
#'
#' @return This function returns a matrix containing a random sample of n
#' vectors drawn form a Tweedie Bayesian network.
#'
#' @examples
#' beta = matrix(0,ncol=5,nrow=5)
#' beta[1,1] = 1;
#' beta[2,1] = 0.5
#' beta[3,1] = -0.5; beta[3,2] = -0.7;beta[3,3] = -1.2
#' beta[4,1] = -2; beta[4,3] = -1
#' beta[5,1] = 0.2; beta[5,4] = 1; beta[5,5] = -2
#' generateSample.tbn(beta, c(1,2,1.5,1,1), c(2,2,3,2,0), 100)
#'
#' @export
generateSample.tbn <- function(beta, phi, p, n) {
  d = length(p)
  X = matrix(0, n, d)
  for (i in 1:d) {
    if (i == 1) {
      PaX = matrix(1, n, 1)
    }
    else{
      PaX = cbind(PaX, X[, i - 1])
    }
    if (p[i] == 0) {
      mu = as.vector(PaX %*% beta[i, 1:i])
    }
    else{
      mu = as.vector(exp(PaX %*% beta[i, 1:i]))
    }

    if (p[i] == 0) {
      X[, i] = rnorm(n, mu, sqrt(phi[i]))
    }
    else if (p[i] == 3) {
      X[, i] = rinvgauss(n, mean = mu, dispersion = phi[i])
    }
    else{
      X[, i] = rtweedie(n,
                        power = p[i],
                        mu = mu,
                        phi = phi[i])
    }
  }
  as.data.frame(X)
}
