# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : decon_hom_err.R
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : Add homoskedastic error argument to 
#
# Input              : None
# Output             : None
# Usage              : 
# 
# References         : None
#
#
# Platform           : R
# Version            : v3.3.0
# Date               : 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Step 1 (of 2): account for homoskedastic error in deconv function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# This function is the same as `deconv` but accounts for 
#   homoskedastic error
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

require('deconvolveR')

deconv_hom_err = function (tau, X, y, Q, P, n = 40, family = c("Poisson", "Normal", 
    "Binomial"), ignoreZero = TRUE, deltaAt = NULL, c0 = 1, scale = TRUE, 
    pDegree = 5, aStart = 1, sd = 1, ...) 
{
    family <- match.arg(family)
    if (missing(Q) && missing(P)) {
        m <- length(tau)
        if (family == "Poisson") {
            if (ignoreZero) {
                supportOfX <- seq_len(n)
                P <- sapply(tau, function(lam) dpois(x = supportOfX, 
                  lambda = lam)/(1 - exp(-lam)))
            }
            else {
                supportOfX <- seq.int(from = 0, to = n - 1)
                P <- sapply(tau, function(w) dpois(x = supportOfX, 
                  lambda = w))
            }
            if (missing(y)) {
                y <- sapply(supportOfX, function(i) sum(X == 
                  i))
            }
            Q <- cbind(1, scale(splines::ns(tau, pDegree), center = TRUE, 
                scale = FALSE))
            Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
        }
        else if (family == "Normal") {
            r <- round(range(X), digits = 1)
            xBin <- seq(from = r[1], to = r[2], length.out = n)
            xBinDropFirst <- xBin[-1]
            xBinDropLast <- xBin[-length(xBin)]
            P <- sapply(tau, function(x) pnorm(q = xBinDropFirst, 
                mean = x, sd = sd) - pnorm(q = xBinDropLast, mean = x, sd = sd))
            intervals <- findInterval(X, vec = xBin)
            y <- sapply(seq_len(n - 1), function(w) sum(intervals == 
                w))
            if (scale) {
                Q1 <- scale(splines::ns(tau, pDegree), center = TRUE, 
                  scale = FALSE)
                Q1 <- apply(Q1, 2, function(w) w/sqrt(sum(w * 
                  w)))
            }
            if (!is.null(deltaAt)) {
                I0 <- as.numeric(abs(tau - deltaAt) < 1e-10)
                Q <- cbind(I0, Q1)
            }
            else {
                Q <- Q1
            }
        }
        else {
            Q <- splines::ns(tau, pDegree)
            if (scale) {
                Q <- scale(Q, center = TRUE, scale = FALSE)
                Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
            }
            P <- sapply(tau, function(w) dbinom(X[, 2], size = X[, 
                1], prob = w))
            y <- 1
        }
    }
    else {
        if (!missing(X) || missing(y) || missing(P) || missing(Q)) {
            stop("P, Q, and y (but not X) must be specified together!")
        }
    }
    p <- ncol(Q)
    pGiven <- length(aStart)
    if (pGiven == 1) {
        aStart <- rep(aStart, p)
    }
    else {
        if (pGiven != p) 
            stop(sprintf("Wrong length (%d) for initial parameter, expecting length (%d)", 
                pGiven, p))
    }
    statsFunction <- function(a) {
        g <- as.vector(exp(Q %*% a))
        g <- g/sum(g)
        G <- cumsum(g)
        f <- as.vector(P %*% g)
        yHat <- if (length(y) == 1 && y == 1) 
            y
        else sum(y) * f
        Pt <- P/f
        W <- g * (t(Pt) - 1)
        qw <- t(Q) %*% W
        ywq <- (yHat * t(W)) %*% Q
        I1 <- qw %*% ywq
        aa <- sqrt(sum(a^2))
        sDot <- c0 * a/aa
        sDotDot <- (c0/aa) * (diag(length(a)) - outer(a, a)/aa^2)
        R <- sum(diag(sDotDot))/sum(diag(I1))
        I2 <- solve(I1 + sDotDot)
        bias <- as.vector(-I2 %*% sDot)
        Cov <- I2 %*% (I1 %*% t(I2))
        Dq <- (diag(g) - outer(g, g)) %*% Q
        bias.g <- Dq %*% bias
        Cov.g <- Dq %*% Cov %*% t(Dq)
        se.g <- diag(Cov.g)^0.5
        D <- diag(length(tau))
        D[lower.tri(D)] <- 1
        Cov.G <- D %*% (Cov.g %*% t(D))
        se.G <- diag(Cov.G)^0.5
        mat <- cbind(tau, g, se.g, G, se.G, bias.g)
        colnames(mat) = c("theta", "g", "SE.g", "G", "SE.G", 
            "Bias.g")
        list(S = R, cov = Cov, cov.g = Cov.g, mat = mat)
    }
    loglik <- function(a) {
        g <- exp(Q %*% a)
        g <- as.vector(g/sum(g))
        f <- as.vector(P %*% g)
        value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
        Pt <- P/f
        W <- g * (t(Pt) - 1)
        qw <- t(Q) %*% W
        aa <- sqrt(sum(a^2))
        sDot <- c0 * a/aa
        if (family == "Binomial") {
            attr(value, "gradient") <- -rowSums(qw) + sDot
        }
        else {
            attr(value, "gradient") <- -(qw %*% y) + sDot
        }
        value
    }
    result <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-10, 
        ...)
    mle <- result$estimate
    stats <- statsFunction(mle)
    list(mle = mle, Q = Q, P = P, S = stats$S, cov = stats$cov, 
        cov.g = stats$cov.g, stats = stats$mat, loglik = loglik, 
        statsFunction = statsFunction)
}

    