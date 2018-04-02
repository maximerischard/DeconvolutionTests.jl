# DeconvolutionTests

[![Build Status](https://travis-ci.org/maximerischard/DeconvolutionTests.jl.svg?branch=master)](https://travis-ci.org/maximerischard/DeconvolutionTests.jl)

[![Coverage Status](https://coveralls.io/repos/maximerischard/DeconvolutionTests.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/maximerischard/DeconvolutionTests.jl?branch=master)

[![codecov.io](http://codecov.io/github/maximerischard/DeconvolutionTests.jl/coverage.svg?branch=master)](http://codecov.io/github/maximerischard/DeconvolutionTests.jl?branch=master)

Test whether two samples come from the same distribution under measurement error. More precisely:

X_i ∼ F_X
ϵ_X_i ∼ Normal(0, σ_X_i)
Xobs = X + ϵ_X

Y_i ∼ F_Y
ϵ_Y_i ∼ Normal(0, σ_Y_i)
Yobs = Y + ϵ

Given the observations Xobs and Yobs, and the measurement errors, σ_X and σ_Y, 
we provide a way to test the null hypothesis that F_X = F_Y.
