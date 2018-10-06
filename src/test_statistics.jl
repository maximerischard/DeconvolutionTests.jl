function numerical_integration(x, y)
    sum(diff(x) .* midpoints(y))
end
function cdf_distance(F_X, F_Y, x, p)
    if p == Inf
        return maximum(abs.(F_X .- F_Y))
    end
    diff_p = abs.(F_X .- F_Y).^p
    return numerical_integration(x, diff_p)^(1/p)
end
function cdf_distance(Xdistr::UnivariateDistribution, Ydistr::UnivariateDistribution, x, p)
    F_X = cdf.(Xdistr, x)
    F_Y = cdf.(Ydistr, x)
    return cdf_distance(F_X, F_Y, x, p)
end
    
function distance_test_statistic(p, method::DeconvolutionMethod, xgrid::AbstractVector)
    function ϕ(X, Y, σ_X, σ_Y)
        ϵ_X = Normal.(0.0, σ_X)
        Xdistr_hat = decon(method, X, ϵ_X)

        ϵ_Y = Normal.(0.0, σ_Y)
        Ydistr_hat = decon(method, Y, ϵ_Y)
        return cdf_distance(Xdistr_hat, Ydistr_hat, xgrid, p)
    end
    return ϕ
end

KS_test_statistic(X, Y, σ_X, σ_Y) = HypothesisTests.ApproximateTwoSampleKSTest(X, Y).δ
