function decon(X::Vector, σ_X::Vector, bw::Float64, grid::Vector; num_t=50, fixup=true)
    ϵ_distr = Normal.(0.0, σ_X)
    Fhat = DeconvolveDistribution.Fhat(grid, X, num_t, bw, ϵ_distr)
    if fixup
        DeconvolveDistribution.fix_CDF!(Fhat)
    end
    return Fhat
end

""" Sample from the CDF """
function cdf_sample(xx, F)
    u = rand() # random uniform
    i = searchsortedfirst(F, u)
    return xx[i]
end
