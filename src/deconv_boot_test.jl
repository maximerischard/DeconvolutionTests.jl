function deconv_boot_test(X::Vector, Y::Vector, σ_X::Vector, σ_Y::Vector, t::Function; niter=100)
    X_and_Y = cat(1, X, Y)
    σ_X_and_Y = cat(1, σ_X, σ_Y)
    t_obs = t(X, Y, σ_X, σ_Y)
    n_X = length(X)
    n_Y = length(Y)
    n_XY = n_X + n_Y
    xx = collect(linspace(minimum(X_and_Y)-1.0, maximum(X_and_Y)+1.0, 1000))
    null_CDF = decon(X_and_Y, σ_X_and_Y, 0.3, xx; fixup=true)
    nabove = 0
    for i in 1:niter
        Xboot = [cdf_sample(xx, null_CDF) for _ in 1:n_X]
        Yboot = [cdf_sample(xx, null_CDF) for _ in 1:n_Y]
        Xtilde = Xboot .+ σ_X.*randn(n_X)
        Ytilde = Yboot .+ σ_Y.*randn(n_Y)
        t_perm = t(Xtilde, Ytilde, σ_X, σ_Y)
        if t_perm > t_obs
            nabove += 1
        end
    end
    return nabove / niter
end
function deconv_boot_stat(X::Vector, Y::Vector, σ_X::Vector, σ_Y::Vector, t::Function; niter=100)
    X_and_Y = cat(1, X, Y)
    σ_X_and_Y = cat(1, σ_X, σ_Y)
    t_obs = t(X, Y, σ_X, σ_Y)
    n_X = length(X)
    n_Y = length(Y)
    n_XY = n_X + n_Y
    t_record = Vector{typeof(t_obs)}(niter)
    xx = collect(linspace(minimum(X_and_Y)-1.0, maximum(X_and_Y)+1.0, 1000))
    null_CDF = decon(X_and_Y, σ_X_and_Y, 0.3, xx; fixup=true)
    nabove = 0
    for i in 1:niter
        Xboot = [cdf_sample(xx, null_CDF) for _ in 1:n_X]
        Yboot = [cdf_sample(xx, null_CDF) for _ in 1:n_Y]
        Xtilde = Xboot .+ σ_X.*randn(n_X)
        Ytilde = Yboot .+ σ_Y.*randn(n_Y)
        t_perm = t(Xtilde, Ytilde, σ_X, σ_Y)
        t_record[i] = t_perm
    end
    return t_record
end
