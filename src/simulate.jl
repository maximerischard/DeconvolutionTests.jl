""" Simulate heteroscedastic data with known errors. """
function sim_data(F_X::Distribution, F_Y::Distribution, 
                  σ_X_distr::Distribution, σ_Y_distr::Distribution, 
                  n_X::Int, n_Y::Int)
    σ_X = rand(σ_X_distr, n_X)
    σ_Y = rand(σ_Y_distr, n_Y)
    X = rand(F_X, n_X)
    Y = rand(F_Y, n_Y)
    ϵX = σ_X .* randn(n_X)
    ϵY = σ_Y .* randn(n_Y)
    Xtilde = X .+ ϵX
    Ytilde = Y .+ ϵY
    return Dict(
        # inputs
        :F_X => F_X,
        :F_Y => F_Y,
        :σ_X_distr => σ_X_distr,
        :σ_Y_distr => σ_Y_distr,
        :n_X => n_X,
        :n_Y => n_Y,
        # outputs
        :σ_X => σ_X,
        :σ_Y => σ_Y,
        :X => X,
        :Y => Y,
        :ϵX => ϵX,
        :ϵY => ϵY,
        :Xtilde => Xtilde,
        :Ytilde => Ytilde,
    )
end

function sim_homo(F_X::Distribution, F_Y::Distribution,
                             σ::Real, n_X::Int, n_Y::Int)
    X = rand(F_X, n_X)
    Y = rand(F_Y, n_Y)
    ϵX = σ .* randn(n_X)
    ϵY = σ .* randn(n_Y)
    Xtilde = X .+ ϵX
    Ytilde = Y .+ ϵY

    σ_X = ones(n_X) .* σ
    σ_Y = ones(n_Y) .* σ

    return Dict(
        # inputs
        :F_X => F_X,
        :F_Y => F_Y,
        :σ => σ,
        :n_X => n_X,
        :n_Y => n_Y,
        # outputs
        :σ_X => σ_X,
        :σ_Y => σ_Y,
        :X => X,
        :Y => Y,
        :ϵX => ϵX,
        :ϵY => ϵY,
        :Xtilde => Xtilde,
        :Ytilde => Ytilde,
    )
end
    
