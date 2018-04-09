library('deconvolveR')
source('../src/decon_hom_err.R')

set.seed(238923) ## for reproducibility


N <- 5000
h = 0.1 			# grid width for tau
theta <- rnorm(N, 0, 4)
X <- theta + rnorm(n = N, 0, 3)
tau <- seq(-20, 20, h)
result <- deconv_hom_err(tau = tau, X = X, family = 'Normal', pDegree = 5, c0 = 0, sd = 3)



plot(result$stats[,1], result$stats[,2]/h, type = 'l', ylim = c(0, 0.1), main = 'Deconvolution with 95% envelopes: error SD == 3',
    xlab = expression(theta), ylab = expression(g(theta)))
lines(result$stats[,1], (result$stats[,2] + 2*result$stats[,3])/h , type = 'l', lty = 2)
lines(result$stats[,1], (result$stats[,2] - 2*result$stats[,3])/h , type = 'l', lty = 2)
lines(tau, dnorm(tau, 0, 4), col = 3)
legend('topright', c('Estimate', 'Truth'), col = c(1, 3), lty = 1)


