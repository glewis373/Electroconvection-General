% This is the new approach to computing the TMatrix.  
% Start with the integral relation
%  q_m(r) = 2 \int_0^\infty k^2 J_m(kr) \int_rin^rout rho J_m(krho) psi2m(rho)
% and change the order of integration 
%  q_m(r) = \int_rin^rout rho psi2m(rho) 2 \int_0^\infty k^2 J_m(kr) J_m(krho) 
%         = \int_rin^rout rho psi2m(rho) Am(r,rho)
% where
%    Am(r,rho) =  2 \int_0^\infty k^2 J_m(kr) J_m(krho) 
%
% There are two steps.  First, we need to approximate Am(r,rho).  
% Because the dk integral
% is over an infinite domain [0,\infty), this is approximated by an integral
% over [0,k_max] with k_max = 50.  (This is the value that Peichun used...)
% The trapezoidal rule is used for the quadrature rule in the k direction.
%
% Once we've approximated Am(r,rho), it's now a matter of approximating
%               \int_rin^rout rho psi2m(rho) Am(r,rho)
% This is done spectrally; construct F = rho.*psi2m.*Am(r,rho) and then
% do the usual (c_int*T_inv*F)/2.
%
% The only problem is: psi2m is embedded in F as part of a pointwise 
% product which is then acted upon by a matrix.  I couldn't figure out % a clever way to express this as qm = T psi2m.  So I did the obvious
% thing of building the matrix T by constructing qm from the standard 
% basis.  Start with qm = <1,0,0,...0> and compute the resutling qm.
% That's the first column of T.  Then take qm = <0,1,0,...0> and compute
% the resulting qm.  That's the second column of T.  And so on.  Not 
% pretty, but it works.
% 
% -------------------------------------------------------------------------
% The function cheb(Nc) in this program is written by Lloyd N. Trefethen.
% -------------------------------------------------------------------------
function T = qm_vs_Psi2m(K,Nc,Nk,k_max, alpha)

[T_inv,c_int,D,DD,x] = ccheb(Nc); % calculating the differentiation matrix for Chebyshev collocation unkowns; x = Chebyshev grids
A = zeros(Nc+1, Nc+1, K+1);

r = x/2 + ( alpha/(1-alpha)+1/2 ); % this is r, with r(1) = r_out and r(Nc+1) = r_in

% Nk = 50000;
% k_max = 50; % These are the values that Peichun used...
dk = k_max/Nk; % Peichun liked dk=.001...
k = [0:dk:k_max]'; % k = [0:dk: k_val()]';
dk = diff(k);
% weights for the integral with respect to k.
delta_k = ([dk;0]+[0;dk])/2;

[rr, kk] = meshgrid(r, k); % Matrix (K, rho') for calculating A_kr = sum( dr' r' J_m(k*r') )

[kk0, rr0] = meshgrid(k, r); % Matrix (rho, K) for calculating q_sumk_rk 
[delta_kk0,rr0] = meshgrid(delta_k,r);

for m = 0:1:K % for m = 0:1:K, for each m
    % A = C*D where
    % C = 2*(kk0.^2).*besselj(m,kk0.*rr0).*delta_kk0;
    % D = besselj(m,kk.*rr);
    A(:,:,m+1) = (2*(kk0.^2).*besselj(m,kk0.*rr0).*delta_kk0)*(besselj(m,kk.*rr));
end

% Now that I've constructed the matrix corresponding to A(rho,r), I need
% to construct T column by column.
I = eye(Nc+1);
for m=0:K
    for i=1:Nc+1
        % test vector is I(:,i).  
        % Operator(I(:,i)) = ith column of T(:,:,m+1).
        T(:,i,m+1) = (c_int*T_inv*(diag(r.*I(:,i))*A(:,:,m+1)))/2;     
    end   
end





