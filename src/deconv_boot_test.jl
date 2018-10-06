function deconv_boot_test(X::Vector, Y::Vector, σ_X::Vector, σ_Y::Vector, t::Function, method::DeconvolutionMethod; niter=100)
    X_and_Y = cat(X, Y, dims=1)
    σ_X_and_Y = cat(σ_X, σ_Y, dims=1)

    ϵ_X = Normal.(0.0, σ_X)
    ϵ_Y = Normal.(0.0, σ_Y)
    ϵ_X_and_Y = Normal.(0.0, σ_X_and_Y)

    t_obs = t(X, Y, σ_X, σ_Y)
    n_X = length(X)
    n_Y = length(Y)
    n_XY = n_X + n_Y
    null_decon = decon(method, X_and_Y, ϵ_X_and_Y)
    nabove = 0
    for i in 1:niter
        Xboot = rand(null_decon, n_X)
        Yboot = rand(null_decon, n_Y)
        Xtilde = Xboot .+ rand.(ϵ_X)
        Ytilde = Yboot .+ rand.(ϵ_Y)
        t_boot = t(Xtilde, Ytilde, σ_X, σ_Y)
        if t_boot > t_obs
            nabove += 1
        end
    end
    return nabove / niter
end

function deconv_boot_stat(X::Vector, Y::Vector, σ_X::Vector, σ_Y::Vector, t::Function, method::DeconvolutionMethod; niter=100)
    X_and_Y = cat(X, Y, dims=1)
    σ_X_and_Y = cat(σ_X, σ_Y, dims=1)

    ϵ_X = Normal.(0.0, σ_X)
    ϵ_Y = Normal.(0.0, σ_Y)
    ϵ_X_and_Y = Normal.(0.0, σ_X_and_Y)

    t_obs = t(X, Y, σ_X, σ_Y)
    n_X = length(X)
    n_Y = length(Y)
    n_XY = n_X + n_Y
    t_record = Vector{typeof(t_obs)}(niter)
    null_decon = decon(method, X_and_Y, ϵ_X_and_Y)
    for i in 1:niter
        Xboot = rand(null_decon, n_X)
        Yboot = rand(null_decon, n_Y)
        Xtilde = Xboot .+ rand.(ϵ_X)
        Ytilde = Yboot .+ rand.(ϵ_Y)
        t_boot = t(Xtilde, Ytilde, σ_X, σ_Y)
        t_record[i] = t_boot
    end
    return t_record
end
