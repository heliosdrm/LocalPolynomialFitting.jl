module LocalPolynomialFitting

export lpparams, lpfit,
       optbandn,
       epanwin

"""Epanechnikov's window"""
epanwin(n::Integer) = 3/4 * (1 - linspace(-1,1,n+2).^2)[2:end-1]

include("lanshammar.jl")

end
