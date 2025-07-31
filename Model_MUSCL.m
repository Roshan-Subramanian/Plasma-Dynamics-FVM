clear all; clc; close all;

%% Parameters
cfl     = 0.5;	% CFL number
tEnd    = 0.15;	% Final time
nx      = 100;  % Number of cells/Elements
n       = 5;	% Number of degrees of freedom in the gas
IC      = 01;	% ~12 Initial value problems are available
limiter ='MM';  % MC, MM, VA.
fluxMth ='ROE'; % LF, ROE, RUS, AUSM, HLLE, HLLC(!).
plot_fig= 1;

% Ratio of specific heats for ideal di-atomic gas
gamma=(n+2)/n;

% Discretize spatial domain
Lx=1; dx=Lx/nx; xc=dx/2:dx:Lx;

% Set IC
[r0,u0,p0] = Euler_IC1d(xc,IC);
E0 = p0./((gamma-1)*r0)+0.5*u0.^2;  % Total Energy
a0 = sqrt(gamma*p0./r0);            % Speed of sound

% Exact solution
[xe,re,ue,pe,ee,te,Me,se] = ...
   EulerExact(r0(1),u0(1),p0(1),r0(nx),u0(nx),p0(nx),tEnd,n);
Ee = pe./((gamma-1)*re)+0.5*ue.^2;

% Set q-array & adjust grid for ghost cells
nx=nx+2; q0=[r0; r0.*u0; r0.*E0]; zero=[0;0;0]; q0=[zero,q0,zero];

% Boundary Conditions in ghost cells
q0(:,1)=q0(:,2); q0(:,nx)=q0(:,nx-1);   % Natural BCs

% Initial time step
lambda0=abs(u0)+a0; dt0=cfl*dx/max(lambda0(:));

% Load IC
q=q0; t=0; it=0; dt=dt0; lambda=lambda0;

% Solver Loop
tic
while t < tEnd
    
    % RK2 1st step
    qs = q - dt*MUSCL_EulerSys(q,max(lambda(:)),gamma,dx,nx,limiter,fluxMth);
    
    qs(:,1)=qs(:,2); qs(:,nx)=qs(:,nx-1);   % Natural BCs
    
    % RK2 2nd step  / update q
    q = (q + qs - dt*MUSCL_EulerSys(qs,max(lambda(:)),gamma,dx,nx,limiter,fluxMth))/2;
    
    q(:,1)=q(:,2); q(:,nx)=q(:,nx-1);   % Natural BCs
        
    % compute flow properties
    r=q(1,:); u=q(2,:)./r; E=q(3,:)./r; p=(gamma-1)*r.*(E-0.5*u.^2); a=sqrt(gamma*p./r); 
    
    % Update dt and time
    lambda=abs(u)+a; dt=cfl*dx/max(lambda(:));
    if t+dt>tEnd; dt=tEnd-t; end
	t=t+dt; it=it+1;
    
    % Plot figure
    if rem(it,10) == 0
        if plot_fig == 1;
            subplot(2,2,1); plot(xc,r(2:nx-1),'.b',xe,re);
            subplot(2,2,2); plot(xc,u(2:nx-1),'.m',xe,ue); 
            subplot(2,2,3); plot(xc,p(2:nx-1),'.k',xe,pe); 
            subplot(2,2,4); plot(xc,E(2:nx-1),'.r',xe,Ee);
            drawnow
        end
    end
end
cputime = toc;

% Remove ghost cells
q=q(:,2:nx-1); nx=nx-2; 

% compute flow properties
r=q(1,:); u=q(2,:)./r; E=q(3,:)./r; p=(gamma-1)*r.*(E-0.5*u.^2);

% Plots results
figure(1);
subplot(2,2,1); plot(xc,r,'ro',xe,re,'-k'); xlabel('x'); ylabel('\rho'); legend('MUSCL','Exact'); 
title('SSP-RK2 TVD-MUSCL Euler Eqns.')
subplot(2,2,2); plot(xc,u,'ro',xe,ue,'-k'); xlabel('x'); ylabel('u'); %legend('MUSCL','Exact');
subplot(2,2,3); plot(xc,p,'ro',xe,pe,'-k'); xlabel('x'); ylabel('p'); %legend('MUSCL','Exact');
subplot(2,2,4); plot(xc,E,'ro',xe,Ee,'-k'); xlabel('x'); ylabel('E'); %legend('MUSCL','Exact');