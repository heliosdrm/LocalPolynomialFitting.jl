# Source: Lanshammar, H. (1982).
# On precision limits for derivatives numerically calculated from noisy data
# Journal of Biomechanics 15(6), 459-470

"""Define parameters for Lanshammar's polynomial filter"""
function lpparams(n, p; window=epanwin, dt::Real=1.)
    U = zeros(2n+1,p+1)
    for k=0:p
      U[:,k+1] = (-n:n).^k
    end
    w = window(2*n+1)
    W = diagm(w./maximum(w))
    V = U'*W*U
    f = map(factorial, 0:p) ./ (-dt).^(0:p)
    Dict(:h => diagm(f)*(V\eye(p+1))*U'*W, :p=>p)
end

"""Apply Lanshammar polynomial filter"""
function lpfit(x::Matrix; h::VecOrMat=[1], p::Integer=0);
y = zeros(size(x)..., p+1)
for k=1:p+1
  y[:,:,k] = filt(collect(h[k,:]),1,x)
end
y
end

function lpfit(x::Vector; h::VecOrMat=[1], p::Integer=0);
y = zeros(length(x), p+1)
for k=1:p+1
  y[:,k] = filt(collect(h[k,:]),1,x)
end
y
end

"""Table of optimal bandwidths"""
optband =
    [NaN  .275  NaN; # p=1
     .53  .275  .33; # p=2
     .53  .600  .33; # p=3
     .82  .600  .70; # p=4
     .82  .920  .70; # p=5
    1.15  .920 1.02; # p=6
    1.15 1.250 1.02; # p=7
    1.48 1.250 1.35; # p=8
    1.48   NaN 1.35] # p=9
    
nmin = 
    [0 1 0; # p=1
     1 1 1; # p=2
     1 2 1; # p=3
     2 2 2; # p=4
     2 3 2; # p=5
     3 3 3; # p=6
     3 4 3; # p=7
     4 4 4; # p=8
     4 0 4] # p=9

function optbandn(p,k,w)
    !(1 <= p <= size(optband)[1]) &&
        error("polynomial order must be between 1 and $(size(optband)[1])")
    !(0 <= k <= size(optband)[2]) &&
        error("the target order of the derivative must be between 0 and $(size(optband)[2])")
    isnan(optband[p,k+1]) && error("value undefined for p==$p, k==$k")
    max(ceil(Integer, optband[p,k+1]/w), nmin[p,k+1])
end
