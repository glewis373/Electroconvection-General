% clear; close all;clc;
global tspan Nc K D2 D riprime dt Re Pr rr
tic

% loading some precomputed input for the time stepper 
% if the parameters (aspect ratio, Fourrier modes and Chebyschev) are
% changed we need to recompute the following output 

load('Tq_vs_psi2_Nc24_K32_rr_0p56_k_decay_value_50.mat')
load('base_state_charge_Nc24_rr_0p56.mat')

Omega = Re*Pr*(1-rr)/rr;


%load('Sim_dat_1','Ra','Ut')
load('Ra621_SL_mv','Ra','Ut')
Ra

% parameter needed for the time stepper of the sheared electroconvection 
tspan = [0 0.1]; %[0 .2];         % time of integration 
K  = 32;                % highest fourier wave
Nc  = 24;               % highest power of the Chebyshev
rr = .56;               % aspect ratio
Re = 0.231; %.249;              % Dimensionless ratio number
Pr = 75.8;              % Dimensionless Prandlt number 
%dt = 1.0e-4;            % time step 
%InitBCs_amp  =  1.0000e-03;

% ==============================================================
% Create Grids and the differentiation Matrices
% ==============================================================
[T_inv,c_int,D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr);

 r_coords = riprime*ones(1,66);
 t_coords = ones(25,1)*linspace(0,2*pi,66);
 xx = r_coords.*cos(t_coords);
 yy = r_coords.*sin(t_coords);


 N_rrr=100; N_ttt=200;  drp= abs(riprime(1)-riprime(end));

 rrr_iprime = (riprime(end):drp/N_rrr:riprime(1));
 ttt_vect=linspace(0,2*pi,N_ttt);
% rrr_coords = rrr_iprime*ones(1,N_ttt);
% ttt_coords = ones(length(rrr_iprime),1)*linspace(0,2*pi,N_ttt);

[rrr_coords,ttt_coords]=meshgrid(rrr_iprime,ttt_vect);

  xxx = rrr_coords.*cos(ttt_coords);
  yyy = rrr_coords.*sin(ttt_coords);

sz_Ut=size(Ut)

cnt=0;

[psi2m,qm,wm,phim,psi2_f,q_f,w_f,phi_f] =  truncvec_2_allmat(Ut(:,1000),Tq_vs_psi2,D2,D,Omega,rr,riprime,Nc,K);

[Phi0,Phi0_r,W0] = base_state_flow(Omega, rr, riprime)

%phi_f = [phi_f phi_f(:,1)];  % phi component
%q_f = [q_f q_f(:,1)];  % q component
%w_f = [w_f w_f(:,1)];  % vorticity component
%psi2_f = [psi2_f psi2_f(:,1)];  % psi2 component

Frames = [1:25:sz_Ut(2)];

% plotting the four physical quantities u = [phi;q;omega;\psi_2]
f1 = figure(1)

for kk = 1:4;

    if (kk == 1) ppp = phi_f + Phi0; tstr='phi';end
    if (kk == 2) ppp = q_f + Q0; tstr='q';end
    if (kk == 3) ppp = w_f + W0; tstr='omega';end
    if (kk == 4) ppp = psi2_f + Psi20; tstr='psi2';end

    ppp = [ppp ppp(:,1)]


    MaxP = max(max(ppp))
    MnP = min(min(ppp))

    nlevs=15;
    levstp=((MaxP-MnP)/nlevs)
    ConLev=MnP:levstp:MaxP


    vvv = interp2(t_coords,r_coords,real(ppp),ttt_coords,rrr_coords,'spline');

%   c1=contourf(xxx,yyy,real(phi_f),ConLev);
    sp = subplot(2,2,kk)
    c1=contourf(xxx,yyy,vvv,ConLev);
         hold on

    plot(xx(1,:),yy(1,:),'k');
    plot(xx(end,:),yy(end,:),'k');

      set(gca,'CLim',[MnP MaxP]);

        hold off

    axis([-2.5 2.5 -2.5 2.5])
    colorbar
    colormap('jet')
   t1=title(tstr,'Interpreter','Latex','FontSize',18);

   pause(0.1)

   set(gcf,'Position',[680 467 553 454]);

end

f1 = figure(2)

for kk = 1:4;

    if (kk == 1) ppp = phi_f; tstr='phi';end
    if (kk == 2) ppp = q_f; tstr='q';end
    if (kk == 3) ppp = w_f; tstr='omega';end
    if (kk == 4) ppp = psi2_f; tstr='psi2';end

    ppp = [ppp ppp(:,1)]


    MaxP = max(max(ppp))
    MnP = min(min(ppp))

    nlevs=15;
    levstp=((MaxP-MnP)/nlevs)
    ConLev=MnP:levstp:MaxP


    vvv = interp2(t_coords,r_coords,real(ppp),ttt_coords,rrr_coords,'spline');

%   c1=contourf(xxx,yyy,real(phi_f),ConLev);
    sp = subplot(2,2,kk)
    c1=contourf(xxx,yyy,vvv,ConLev);
         hold on

    plot(xx(1,:),yy(1,:),'k');
    plot(xx(end,:),yy(end,:),'k');

      set(gca,'CLim',[MnP MaxP]);

        hold off

    axis([-2.5 2.5 -2.5 2.5])
    colorbar
    colormap('jet')
   t1=title(tstr,'Interpreter','Latex','FontSize',18);

   pause(0.1)

   set(gcf,'Position',[680 467 553 454]);

end

