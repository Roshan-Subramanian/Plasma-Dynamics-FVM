function [Q_values] = MUSCL_Scheme(Q,gamma,delx,N)

% Output
Q_values = zeros(3,N);
dQ = zeros(3,N); 
% Flux
F=zeros(3,N-1); 
% Left step in MUSCL
Q_L=zeros(3,N-1); 
% Right step in MUSCL
Q_R=zeros(3,N-1);

for i=1:3
    for j=2:N-1
        % Minmod limiter
        % Slope
        dQ_R = (Q(i,j+1)-Q(i,j))/delx;
        dQ_L = (Q(i,j) - Q(i,j-1))/delx;
        dQ(i,j) = minmod(dQ_L,dQ_R);
    end
end

% Linear Extrapolation till j-1/2, and j+1/2
% 2nd order, and fully upwind
for j=2:N-2
    Q_L(:,j) = Q(:,j) + dQ(:,j)*delx/2;
    Q_R(:,j) = Q(:,j+1) - dQ(:,j+1)*delx/2;
end

%% Rieman Solver to determine edge flux through cell interface

% Flux for the inner cells
for j=2:N-2
    F(:,j) = ROEflux(Q_L(:,j),Q_R(:,j),gamma);
    Q_values(:, j ) = Q_values(:, j ) + F(:,j)/delx;
    Q_values(:,j+1) = Q_values(:,j+1) - F(:,j)/delx;
end
 
% Flux for the left wall
Q_R(:,1)=Q(:,2)-dQ(:,2)*delx/2;    
Q_L(:,1)=Q_R(:,1);
F(:,1) = ROEflux(Q_L(:,1),Q_R(:,1),gamma);
Q_values(:,2) = Q_values(:,2) - F(:,1)/delx;

% Flux for the right wall
Q_L(:,N-1)=Q(:,N-1)+dQ(:,N-1)*delx/2;      
Q_R(:,N-1) = Q_L(:,N-1);
F(:,N-1) = ROEflux(Q_L(:,N-1),Q_R(:,N-1),gamma);
Q_values(:,N-1) = Q_values(:,N-1) + F(:,N-1)/delx;

end

 function theta_minmod = minmod(dQ_L,dQ_R)
    if sign(dQ_L) ~= sign(dQ_R) || sign(dQ_L)==0 && sign(dQ_R)==0
        theta_minmod = 0;
    else
        min_dQ = min(abs(dQ_L),abs(dQ_R));
        if min_dQ == abs(dQ_L)
            signed = sign(dQ_L);
            theta_minmod = signed*min_dQ;
        else
            signed = sign(dQ_R);
            theta_minmod = signed*min_dQ;
        end
     end
end

% ROE flux
function Flux_values = ROEflux(Q_L,Q_R,gamma)
    
    % Left state
    rho_L = Q_L(1);
    u_L = Q_L(2)./rho_L;
    E_L = Q_L(3)./rho_L;
    p_L = (gamma-1)*( Q_L(3) - rho_L*u_L*u_L/2 );
    a_L = sqrt(gamma*p_L/rho_L);
    H_L = ( Q_L(3) + p_L )./rho_L;
    
    % Right state
    rho_R = Q_R(1);
    u_R = Q_R(2)./rho_R;
    E_R = Q_R(3)./rho_R;
    p_R = (gamma-1)*( Q_R(3) - rho_R*u_R*u_R/2 );
    a_R = sqrt(gamma*p_R/rho_R);
    H_R = ( Q_R(3) + p_R )./rho_R;
    
    % First compute the Roe Averages
    RT = sqrt(rho_R/rho_L);
    rho = RT*rho_L;
    u = (u_L+RT*u_R)/(1+RT);
    H = (H_L+RT*H_R)/(1+RT);
    a = sqrt((gamma-1)*(H-u*u/2) );
    
    % Differences in primitive variables.
    dr = rho_R - rho_L;
    du = u_R - u_L;
    dP = p_R - p_L;
    
    % Wave strength (Characteristic Variables).
    CV = [(dP-rho*a*du)/(2*a^2); -(dP/(a^2)-dr); (dP+rho*a*du)/(2*a^2)];
    
    % Absolute values of the wave speeds (Eigenvalues)
    EV = [ abs(u-a); abs( u ); abs(u+a) ];


    % Right eigenvectors
    R = [  1  ,  1  ,  1 ;u-a ,  u  , u+a ;H-u*a,u^2/2,H+u*a];
   
    % Compute the average flux.
    F_L=[rho_L.*u_L; rho_L.*u_L.^2+p_L; u_L.*(rho_L.*E_L+p_L)];
    F_R=[rho_R.*u_R; rho_R.*u_R.^2+p_R; u_R.*(rho_R.*E_R+p_R)];

    % Add the matrix dissipation term to complete the Roe flux.
    Flux_values = ( F_L + F_R  - R*(EV.*CV))/2;

end
