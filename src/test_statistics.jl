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
function distance_test_statistic(p, bw, xgrid)
    function ϕ(X, Y, σ_X, σ_Y)
        F_X_hat = decon(X, σ_X, bw, xgrid)
        F_Y_hat = decon(Y, σ_Y, bw, xgrid)
        return cdf_distance(F_X_hat, F_Y_hat, xgrid, p)
    end
    return ϕ
end
