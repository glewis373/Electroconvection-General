function [psi2m,qm,wm,phim,psi2,q,w,phi] =  truncvec_2_allmat(U_t,Tq_vs_psi2,D2,D,Omega,rr,riprime,Nc,K)
len = (Nc+1)*(2*K+1);

psi2 = reshape(U_t(1:len),(Nc+1),(2*K+1));
phi    = reshape(U_t(len+1:2*len),(Nc+1),(2*K+1));


[Phi0,Phi0_r,W0] = base_state_flow(Omega, rr, riprime);
Phi0 = Phi0*ones(1,2*K+1);
%Utheta0 = -Phi0_r*ones(1,2*K+1);
W0 = W0*ones(1,2*K+1);
%Flux0 = -Phi0_r;
%
Phi = Phi0 + phi;

% to get initial data for charge density
psi2m = my_fft(psi2,2*K+1);
for m=0:1:K              %  Assign \hat{q}_m in accordance with the values of  \hat{psi2}_m.
    qm(:, m+1)=Tq_vs_psi2(:,:, m+1)*psi2m(:, m+1);
end
q = my_ifft(qm,2*K+1);


Phim = my_fft(Phi,2*K+1);
phim = my_fft(phi,2*K+1);

for m=0:1:K
    % this is the Laplacian
    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
   wm(:,m+1) = -Lop*phim(:,m+1);
end
w = my_ifft(wm,2*K+1);



end


%for m=0:1:K
%
%
%if m == 0
%    % have boundary conditions on the zero mode
%    % here is the (negative) of the velocity at the inner radius
%    h1 = - in_speed;
%    if phi_BC == 1
%        % here is the total flux
%        % g1 = Chebyshev_Int(Nc,flux_new)/2;
%        g1 = dot(c_int,T_inv*flux_new)/2;
%    else
%        g1 = 0;
%    end
%    % because the boundary conditions are imposed in Fourier space, not in
%    % real space, we need to multiply by the factor the the FFT brings to
%    % the party.
%    h1 = (2*K+1)*h1;
%    g1 = (2*K+1)*g1;
%else
%    h1 = 0;
%    g1 = 0;
%end
%g2 = 0;
%h2 = 0;
%
%BCs4w (:, m+1) = Mi(:,:,m+1)\( [h1;h2] - [2*D(Nc+1, :)*Phim2(:,m+1); 2*D(1, :)*Phim2(:, m+1)] );
%% define the solution as the sum of the step 1 and step 2 problems
%%Wm_new(:,m+1)=Wm2(:,m+1)+BCs4w(1,m+1)*wm1(:,m+1)+BCs4w(2,m+1)*wm3(:,m+1);
%Phim2(:,m+1)=Phim_in(:,m+1)-BCs4w(1,m+1)*phim1(:,m+1)-BCs4w(2,m+1)*phim3(:,m+1);
%
%
%    % this is the Laplacian
%    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
%
%%% Now compute the vorticity from Phim...
%A_full = Lop;
%A_trunc = A_full(2:Nc,2:Nc);
%% here are the boundary conditions on the potential
%%Phim2(1,m+1)=g2;
%%Phim2(Nc+1,m+1)=g1;
%%B = -Wm2(2:Nc,m+1);
%%RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[Phim2(1,m+1);Phim2(Nc+1,m+1)];
%%Phim2(2:Nc,m+1) = A_trunc\RHS;
%
%% compute the coefficients needed to solve the "step 2 problem"
%BCs4w (:, m+1) = Mi(:,:,m+1)\( [h1;h2] - [2*D(Nc+1, :)*Phim2(:,m+1); 2*D(1, :)*Phim2(:, m+1)] );
%% define the solution as the sum of the step 1 and step 2 problems
%%Wm_new(:,m+1)=Wm2(:,m+1)+BCs4w(1,m+1)*wm1(:,m+1)+BCs4w(2,m+1)*wm3(:,m+1);
%Phim2(:,m+1)=Phim_in(:,m+1)-BCs4w(1,m+1)*phim1(:,m+1)-BCs4w(2,m+1)*phim3(:,m+1);
%
%
%
%   wm(:,m+1) = -Lop*phim(:,m+1);
%%end
%w = my_ifft(wm,2*K+1);
%
%
%
%
%%W = my_ifft(Wm,2*K+1);
%%w = W - W0;
%%wm = my_fft(w,2*K+1);
%
%%% Now compute the potential from wm2...
%%A_full = Lop;
%%A_trunc = A_full(2:Nc,2:Nc);
%%% here are the boundary conditions on the potential
%%Phim2(1,m+1)=g2;
%%Phim2(Nc+1,m+1)=g1;
%%B = -Wm2(2:Nc,m+1);
%%RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[Phim2(1,m+1);Phim2(Nc+1,m+1)];
%%Phim2(2:Nc,m+1) = A_trunc\RHS;
%
%
%%
