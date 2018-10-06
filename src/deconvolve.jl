# function decon(method::DeconvolutionMethod, X::Vector, σ_X::Vector)
    # ϵ_distr = Normal.(0.0, σ_X)
    # deconvolved = decon(method, X, ϵ_distr)
    # return deconvolved
# end

# """ Sample from the CDF """
# function cdf_sample(xx, F)
    # u = rand() # random uniform
    # i = searchsortedfirst(F, u)
    # return xx[i]
# end
