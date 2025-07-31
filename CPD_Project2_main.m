
%% Project-2

%% Instructions to run the code

%This shows the Euler simulation
%The files needed are CPD_Project2_main.m, and MUSCL_Scheme.m.

%You will get the results directly by just running the CPD_Project2_main.m file.

%% MUSCL scheme Finite Volume model


clear all
clc;

% Initial Conditions from SOD paper
% p1>p2, rho1>rho2, u1=u2=0



%      Region 1             Region 2
%       p1
%-------------------|
%       rho1        |
%-------------------|         p2
%                   |-----------------------
%                   |      rho2
%                   |-----------------------
%                   |       u2=0
%     u1=0          |       
%-------------------------------------------



% Region 1
rho1 = 1.0;
p1 = 1.0;
u1 = 0.0;

% Region 2
rho2 = 0.125;
p2 = 0.1;
u2 = 0.0;

rho = [rho1 rho2];
p = [p1 p2];
u = [u1 u2];

gamma = 1.4;  % Ratio of specific heat
n_cells = 102; % Number of cells
cfl = 0.9; % cfl number
t_end = 0.2; % End time

% Discretizing spatial domain
delx = 0.01; % From SOD paper
Lx = 1; % Length
x = 0.001:delx:1;

% Splitting the regions
rho0 = zeros(size(x)); 
u0 = zeros(size(x)); 
p0 = zeros(size(x));

% Parameters of regions dimensions
x_middle = (x(end)-x(1))/2;
Left = find(x<=x_middle);
Right = find(x>x_middle);

rho0(Left) = rho(1); % region 1
rho0(Right) = rho(2); % region 2

u0(Left) = u(1); % region 1
u0(Right) = u(2); % region 2

p0(Left) = p(1); % region 1
p0(Right) = p(2); % region 2

% Total Energy
%E = (p0./((gamma-1)))+(0.5*((rho0.*u0).^2)./(rho0));
E = p0./((gamma-1)*rho0)+0.5*u0.^2;
% Speed of sound
a = sqrt(gamma*p0./rho0);
shock_speed = abs(u0)+a;


% Time step
delt0 = cfl*delx/max(shock_speed(:));

% Euler Equation - Q array
Q = [rho0; rho0.*u0; rho0.*E];
zero_array = [0;0;0];
Q = [zero_array,Q,zero_array];

% Boundary Conditions
Q(:,1)=Q(:,2); 
Q(:,n_cells)=Q(:,n_cells-1); 

% Time IC
t = 0; delt = delt0;

while t<t_end
      
    Q1 = Q-delt*MUSCL_Scheme(Q,gamma,delx,n_cells); 
    Q1(:,1) = Q1(:,2);
    Q1(:,n_cells) = Q1(:,n_cells-1);
    
    Q = Q+Q1-delt*MUSCL_Scheme(Q1,gamma,delx,n_cells); 
    Q(:,1)=Q(:,2);
    Q(:,n_cells)=Q(:,n_cells-1);
    
    rho = Q(1,:); u=Q(2,:)./rho; E=Q(3,:)./rho;
    %p=(gamma-1)*rho.*(E-(0.5*((rho.*u).^2)./(rho))); 
    p=(gamma-1)*rho.*(E-0.5*u.^2);
    a=sqrt(gamma*p./rho);
    
    shock_speed = abs(u)+a; 
    % dynamic time stepping using cfl
    delt = cfl*delx/max(shock_speed(:));
     if t+delt>t_end; delt=t_end-t; end
        t=t+delt;
    
end

% Neglecting the left, and right walls
Q=Q(:,2:n_cells-1); n_cells=n_cells-2; 

% flow properties
rho=Q(1,:); u=Q(2,:)./rho; E=Q(3,:)./rho; 
rho = rho/max(Q(1,:)); 
p=(gamma-1)*rho.*(E-0.5*u.^2);

% Plots results
figure(1);
subplot(2,2,1); 
plot(x,rho,'r','Linewidth',2); 
xlabel('x'); ylabel('\rho');
title('Density plot - Euler Equation');
subplot(2,2,2);
plot(x,u,'b','Linewidth',2); 
xlabel('x'); ylabel('u');
title('Velocity plot - Euler Equation');
subplot(2,2,3);
plot(x,p,'c','Linewidth',2); 
xlabel('x'); ylabel('P');
title('Pressure plot - Euler Equation');
subplot(2,2,4);
plot(x,E,'g','Linewidth',2); 
xlabel('x'); ylabel('E');
title('Energy plot - Euler Equation');






