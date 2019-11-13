
function P = updatePP(Y,mu,A, X)
G = X';
Q = (A - Y/mu)';
W = G'* Q + eps;

[U,S,V] = svd (W,'econ'); 

PT = U*V';
P = PT';
end
