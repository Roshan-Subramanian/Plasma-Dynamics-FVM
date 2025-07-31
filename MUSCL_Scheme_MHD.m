function [Q_values] = MUSCL_Scheme_MHD(Q,gamma,delx,N)

% Output
Q_values = zeros(7,N);
dQ = zeros(7,N); 
% Flux
F=zeros(7,N-1); 
% Left step in MUSCL
Q_L=zeros(7,N-1); 
% Right step in MUSCL
Q_R=zeros(7,N-1);

for i=1:7
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
    %F(:,j) = ROEflux(Q_L(:,j),Q_R(:,j),gamma);
    F(:,j) = HLLflux(Q_L(:,j),Q_R(:,j),gamma);
    Q_values(:, j ) = Q_values(:, j ) + F(:,j)/delx;
    Q_values(:,j+1) = Q_values(:,j+1) - F(:,j)/delx;
end
 
% Flux for the left wall
Q_R(:,1)=Q(:,2)-dQ(:,2)*delx/2;    
Q_L(:,1)=Q_R(:,1);
%F(:,1) = ROEflux(Q_L(:,1),Q_R(:,1),gamma);
F(:,1) = HLLflux(Q_L(:,1),Q_R(:,1),gamma);
Q_values(:,2) = Q_values(:,2) - F(:,1)/delx;

% Flux for the right wall
Q_L(:,N-1)=Q(:,N-1)+dQ(:,N-1)*delx/2;      
Q_R(:,N-1) = Q_L(:,N-1);
%F(:,N-1) = ROEflux(Q_L(:,N-1),Q_R(:,N-1),gamma);
F(:,1) = HLLflux(Q_L(:,1),Q_R(:,1),gamma);
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
    % Permeability of free space
    mu0 = 4*pi*1e-07;
    Bx = 0.75;
    % Left state
    rho_L = Q_L(1);
    vx_L = Q_L(2)./rho_L;
    vy_L = Q_L(3)./rho_L;
    vz_L = Q_L(4)./rho_L;
    %Bx_L = Q_L(5)./rho_L;
    By_L = Q_L(5)./rho_L;
    Bz_L = Q_L(6)./rho_L;
    E_L = Q_L(7)./rho_L;
    % Total velocity
    vt_L = sqrt((vx_L).^2+(vy_L).^2+(vz_L).^2);
    % Total Magnetic field
    Bt_L = sqrt((Bx).^2+(By_L).^2+(Bz_L).^2);
    B_perp_L = sqrt((By_L).^2+(Bz_L).^2);
    %p_L = (gamma-1)*(Q_L(7) - rho_L*vt_L*vt_L/2 - Bt_L/(2*mu0));
    %p_L= 2/3.*(E_L-0.5.*rho_L.*vt_L.*vt_L-Bt_L.^2/(2.*mu0));
    p_L = rho_L.*E_L+(Bt_L.^2/(2.*mu0));
    a_L = sqrt(gamma*p_L/rho_L);
    H_L = ( Q_L(7) + p_L )./rho_L;

    
    
    % Right state
    rho_R = Q_R(1);
    vx_R = Q_R(2)./rho_R;
    vy_R = Q_R(3)./rho_R;
    vz_R = Q_R(4)./rho_R;
    %Bx_R = Q_R(5)./rho_R;
    By_R = Q_R(5)./rho_R;
    Bz_R = Q_R(6)./rho_R;
    E_R = Q_R(7)./rho_R;
    % Total velocity
    vt_R = sqrt((vx_R).^2+(vy_R).^2+(vz_R).^2);
    % Total Magnetic field
    Bt_R = sqrt((Bx).^2+(By_R).^2+(Bz_R).^2);
    B_perp_R = sqrt((By_R).^2+(Bz_R).^2);
    %p_R = (gamma-1)*(Q_R(7) - rho_R*vt_R*vt_R/2 - Bt_R/(2*mu0));
    %p_R = 2/3.*(E_R-0.5.*rho_R.*vt_R.*vt_R-Bt_R.^2/(2.*mu0));
    p_R = rho_R.*E_R+(Bt_R.^2/(2.*mu0));
    a_R = sqrt(gamma*p_R/rho_R);
    H_R = ( Q_R(7) + p_R )./rho_R;
    
    % Roe Averages
    rho_bar = sqrt(rho_L).*sqrt(rho_R);
    %rho_bar = ((sqrt(rho_L).*rho_L)+(sqrt(rho_R).*rho_R))/sqrt(rho_L)+sqrt(rho_R);
    u_bar = ((sqrt(rho_L).*vx_L)+(sqrt(rho_R).*vx_R))/sqrt(rho_L)+sqrt(rho_R);
    v_bar = ((sqrt(rho_L).*vy_L)+(sqrt(rho_R).*vy_R))/sqrt(rho_L)+sqrt(rho_R);
    w_bar = ((sqrt(rho_L).*vz_L)+(sqrt(rho_R).*vz_R))/sqrt(rho_L)+sqrt(rho_R);
    vt_bar = ((sqrt(rho_L).*vt_L)+(sqrt(rho_R).*vt_R))/sqrt(rho_L)+sqrt(rho_R);
    B_bar_yz = ((sqrt(rho_R).*B_perp_L)+(sqrt(rho_L).*B_perp_R))/sqrt(rho_L)+sqrt(rho_R);
    B_bar = ((sqrt(rho_R).*Bt_L)+(sqrt(rho_L).*Bt_R))/sqrt(rho_L)+sqrt(rho_R);
    By_bar = ((sqrt(rho_R).*By_L)+(sqrt(rho_L).*By_R))/sqrt(rho_L)+sqrt(rho_R);
    Bz_bar = ((sqrt(rho_R).*Bz_L)+(sqrt(rho_L).*Bz_R))/sqrt(rho_L)+sqrt(rho_R);
    H_bar = ((sqrt(rho_L).*H_L)+(sqrt(rho_R).*H_R))/sqrt(rho_L)+sqrt(rho_R);
    a_bar = ((sqrt(rho_L).*a_L)+(sqrt(rho_R).*a_R))/sqrt(rho_L)+sqrt(rho_R);
    
    
%     del21 = -u_bar.^2+(2-gamma).*X+((gamma-1)./2).*vt_bar.^2;
%     del22 = 2.*u_bar-(gamma-1).*u_bar;
%     del23 = -(gamma-1).*v_bar;
%     del24 = -(gamma-1).*w_bar;
%     del25 = (2-gamma).*By_bar;
%     del26 = (2-gamma).*Bz_bar;
%     del27 = gamma-1;
%     
%     del71 = -u_bar.*H_bar+u_bar.*(del21+u_bar.^2)+(Bx./rho_bar).*(B_bar);
%     del72 = H_bar+u_bar.*(del22-2.*u_bar)-(Bx.^2/rho_bar);
%     del73 = u_bar.*del23-(Bx./rho_bar).*By_bar;
%     del74 = u_bar.*del24-(Bx./rho_bar).*Bz_bar;
%     del75 = u_bar.*del25-Bx.*v_bar;
%     del76 = u_bar.*del26-Bx.*w_bar;
%     del77 = u_bar+u_bar.*del27;
%     
    % Roe matrix
%     Roe_Matrix = [0,1,0,0,0,0,0;del21,del22,del23,del24,del25,del26,del27;...
%         -u_bar.*v_bar,v_bar,u_bar,0,-Bx,0,0;-u_bar.*w_bar,w_bar,0,u_bar,0,-Bx,0;...
%         ((-By_bar./rho_bar).*u_bar)+((Bx./rho_bar).*v_bar),(By_bar./rho_bar),-Bx./rho_bar,0,u_bar,0,0;...
%         ((-Bz_bar./rho_bar).*u_bar)+((Bx./rho_bar).*w_bar),Bz_bar./rho_bar,0,-Bx./rho_bar,0,u_bar,0;...
%         del71,del72,del73,del74,del75,del76,del77];
    
    
    %Alfven speed
    ca_L = sqrt(Bt_L.^2/(4*pi*rho_L));
    ca_R = sqrt(Bt_R.^2/(4*pi*rho_R));
    cax_L = sqrt(Bx.^2/(4*pi*rho_L));
    cax_R = sqrt(Bx.^2/(4*pi*rho_R));
    
    % fast magnetosonic speed
    cf_L = sqrt(1/2.*(a_L.^2+ca_L.^2)+1/2.*sqrt((a_L.^2+ca_L.^2).^2-4*a_L.^2*cax_L.^2));
    cf_R = sqrt(1/2.*(a_R.^2+ca_R.^2)+1/2.*sqrt((a_R.^2+ca_R.^2).^2-4*a_R.^2*cax_R.^2));
    
    % slow magnetosonic speed
    cs_L = sqrt(1/2.*(a_L.^2+ca_L.^2)-1/2.*sqrt((a_L.^2+ca_L.^2).^2-4*a_L.^2*cax_L.^2));
    cs_R = sqrt(1/2.*(a_R.^2+ca_R.^2)-1/2.*sqrt((a_R.^2+ca_R.^2).^2-4*a_R.^2*cax_R.^2));
    
    cf_bar = ((sqrt(rho_L).*cf_L)+(sqrt(rho_R).*cf_R))/sqrt(rho_L)+sqrt(rho_R);
    cs_bar = ((sqrt(rho_L).*cs_L)+(sqrt(rho_R).*cs_R))/sqrt(rho_L)+sqrt(rho_R);
    ca_bar = ((sqrt(rho_L).*ca_L)+(sqrt(rho_R).*ca_R))/sqrt(rho_L)+sqrt(rho_R);
    
    B_perp_bar = sqrt(By_bar.^2+Bz_bar.^2);
    %Beta_yz = B_bar_yz./B_perp_bar;
    Beta_y = By_bar./B_perp_bar;
    Beta_z = Bz_bar./B_perp_bar;
    S = sign(Bx);
    alpha_f = sqrt((a_bar.^2-cs_bar.^2)/(cf_bar.^2-cs_bar.^2));
    alpha_s = sqrt((cf_bar.^2-a_bar.^2)/(cf_bar.^2-cs_bar.^2));

    % Differences in primitive variables.
    Del_Bt = Bt_R-Bt_L;
    Del_By = By_R-By_L;
    Del_Bz = Bz_R-Bz_L;
    Del_rho = rho_R-rho_L;
    Del_p = p_R-p_L;
    Del_u = vx_R-vx_L;
    Del_v = vy_R-vy_L;
    Del_w = vz_R-vz_L;
    
    X = (Del_Bt).^2/(2*(sqrt(rho_L)+sqrt(rho_R)).^2);
    
     % Wave strength (Characteristic Variables)
    CV_u_bar = (a_bar.^2-X).*Del_rho-Del_p;
    CV_u_bar_plus_ca_bar = 1/2*(-Beta_y.*Del_w+Beta_z.*Del_v+(S./sqrt(rho_bar))...
        .*(Beta_y.*Del_Bz-Beta_z.*Del_By));
    CV_u_bar_minus_ca_bar = 1/2*(Beta_y.*Del_w-Beta_z.*Del_v+(S./sqrt(rho_bar))...
        .*Beta_y.*Del_Bz-Beta_z.*Del_By);
    CV_u_bar_plus_cs_bar = 1/2*(alpha_s.*(X.*Del_rho+Del_p)+rho_bar.*alpha_f.*cf_bar.*S...
        .*(Beta_y.*Del_v+Beta_z.*Del_w)+rho_bar.*alpha_s.*cs_bar.*Del_u-sqrt(rho_bar)...
        .*alpha_f.*a_bar.*(Beta_y.*Del_By+Beta_z.*Del_Bz));
    CV_u_bar_minus_cs_bar = 1/2*(alpha_s.*(X.*Del_rho+Del_p)-rho_bar.*alpha_f.*cf_bar.*S...
        .*(Beta_y.*Del_v+Beta_z.*Del_w)-rho_bar.*alpha_s.*cs_bar.*Del_u-sqrt(rho_bar)...
        .*alpha_f.*a_bar.*(Beta_y.*Del_By+Beta_z.*Del_Bz));
    CV_u_bar_plus_cf_bar = 1/2*(alpha_f.*(X.*Del_rho+Del_p)-rho_bar.*alpha_s.*cs_bar.*S...
        .*(Beta_y.*Del_v+Beta_z.*Del_w)+rho_bar.*alpha_f.*cf_bar.*Del_u+sqrt(rho_bar)...
        .*alpha_s.*a_bar.*(Beta_y.*Del_By+Beta_z.*Del_Bz));
    CV_u_bar_minus_cf_bar = 1/2*(alpha_f.*(X.*Del_rho+Del_p)+rho_bar.*alpha_s.*cs_bar.*S...
        .*(Beta_y.*Del_v+Beta_z.*Del_w)+rho_bar.*alpha_f.*cf_bar.*Del_u+sqrt(rho_bar)...
        .*alpha_s.*a_bar.*(Beta_y.*Del_By+Beta_z.*Del_Bz));
    
    CV = [CV_u_bar_minus_cf_bar;CV_u_bar_minus_ca_bar;CV_u_bar_minus_cs_bar;...
        CV_u_bar;CV_u_bar_plus_cs_bar;CV_u_bar_plus_ca_bar;CV_u_bar_plus_cf_bar];
    
    % Absolute values of the wave speeds (Eigenvalues)
    ev = [ abs(u_bar-cf_bar);abs(u_bar-ca_bar);abs(u_bar-cs_bar);abs(u_bar);...
        abs(u_bar+cs_bar);abs(u_bar+ca_bar);abs(u_bar+cf_bar)]; 

    % Right eigenvectors
    R_u_bar = 1/a_bar.^2.*[1,u_bar,v_bar,w_bar,0,0,(vt_bar.^2/2)+((gamma-2)/(gamma-1)).*X];
    R_u_bar_plus_ca_bar = [0,0,rho_bar.*Beta_z,-rho_bar.*Beta_y,-S.*sqrt(rho_bar).*Beta_z,...
        S.*sqrt(rho_bar).*Beta_y,rho_bar.*(v_bar.*Beta_z-w_bar.*Beta_y)];
    R_u_bar_minus_ca_bar = [0,0,-rho_bar.*Beta_z,rho_bar.*Beta_y,-S.*sqrt(rho_bar).*Beta_z,...
        S.*sqrt(rho_bar).*Beta_y,-rho_bar.*(v_bar.*Beta_z-w_bar.*Beta_y)];
    R_u_bar_plus_cs_bar = 1./(rho_bar.*a_bar.^2).*[rho_bar.*alpha_s,rho_bar.*alpha_s.*(u_bar+cs_bar),...
        rho_bar.*(alpha_s.*v_bar+alpha_f.*cf_bar.*Beta_y.*S),rho_bar.*(alpha_s.*w_bar+alpha_f.*cf_bar.*Beta_z.*S),...
        -sqrt(rho_bar).*alpha_f.*a_bar.*Beta_y,-sqrt(rho_bar).*alpha_f.*a_bar.*Beta_z,...
        rho_bar.*alpha_s.*(H_bar-(B_bar.^2./rho_bar)+u_bar.*cs_bar)+rho_bar.*alpha_f.*cf_bar.*S.*(v_bar.*Beta_y...
        +w_bar.*Beta_z)-sqrt(rho_bar).*alpha_f.*a_bar.*B_perp_bar];
    R_u_bar_minus_cs_bar = 1./(rho_bar.*a_bar.^2).*[rho_bar.*alpha_s,rho_bar.*alpha_s.*(u_bar-cs_bar),...
        rho_bar.*(alpha_s.*v_bar-alpha_f.*cf_bar.*Beta_y.*S),rho_bar.*(alpha_s.*w_bar-alpha_f.*cf_bar.*Beta_z.*S),...
        -sqrt(rho_bar).*alpha_f.*a_bar.*Beta_y,-sqrt(rho_bar).*alpha_f.*a_bar.*Beta_z,...
        rho_bar.*alpha_s.*(H_bar-(B_bar.^2./rho_bar)-u_bar.*cs_bar)-rho_bar.*alpha_f.*cf_bar.*S.*(v_bar.*Beta_y...
        +w_bar.*Beta_z)-sqrt(rho_bar).*alpha_f.*a_bar.*B_perp_bar];
    R_u_bar_plus_cf_bar = 1./(rho_bar.*a_bar.^2).*[rho_bar.*alpha_f,rho_bar.*alpha_f.*(u_bar+cf_bar),...
        rho_bar.*(alpha_f.*v_bar-alpha_s.*cs_bar.*Beta_y.*S),rho_bar.*(alpha_f.*w_bar-alpha_s.*cs_bar.*Beta_z.*S),...
        sqrt(rho_bar).*alpha_s.*a_bar.*Beta_y,sqrt(rho_bar).*alpha_s.*a_bar.*Beta_z,...
        rho_bar.*alpha_f.*(H_bar-(B_bar.^2./rho_bar)+u_bar.*cf_bar)-rho_bar.*alpha_s.*cs_bar.*S.*(v_bar.*Beta_y...
        +w_bar.*Beta_z)-sqrt(rho_bar).*alpha_s.*a_bar.*B_perp_bar];
    R_u_bar_minus_cf_bar = 1./(rho_bar.*a_bar.^2).*[rho_bar.*alpha_f,rho_bar.*alpha_f.*(u_bar-cf_bar),...
        rho_bar.*(alpha_f.*v_bar+alpha_s.*cs_bar.*Beta_y.*S),rho_bar.*(alpha_f.*w_bar+alpha_s.*cs_bar.*Beta_z.*S),...
        sqrt(rho_bar).*alpha_s.*a_bar.*Beta_y,sqrt(rho_bar).*alpha_s.*a_bar.*Beta_z,...
        rho_bar.*alpha_f.*(H_bar-(B_bar.^2./rho_bar)-u_bar.*cf_bar)+rho_bar.*alpha_s.*cs_bar.*S.*(v_bar.*Beta_y...
        +w_bar.*Beta_z)-sqrt(rho_bar).*alpha_s.*a_bar.*B_perp_bar];
   
    
   
    R = [R_u_bar_minus_cf_bar;R_u_bar_minus_ca_bar;R_u_bar_minus_cs_bar;R_u_bar;
        R_u_bar_plus_cs_bar;R_u_bar_plus_ca_bar;R_u_bar_plus_cf_bar];
    
   
    F_L=[rho_L.*vt_L;rho_L.*vt_L.^2+(Bx.^2./mu0)+p_L+(Bt_L.^2/(2.*mu0));...
         rho_L.*vx_L.*vy_L-(Bx.*By_L)./mu0;rho_L.*vx_L.*vz_L-(Bx.*Bz_L)./mu0;...
         vx_L.*By_L-Bx.*vy_L;vx_L.*Bz_L-Bx.*vz_L;...
         (E_L+p_L+(Bt_L.^2./(2.*mu0))).*vx_L-((Bt_L.*vt_L)./mu0).*Bx];
        
    F_R=[rho_R.*vt_R;rho_R.*vt_R.^2+(Bx.^2./mu0)+p_R+(Bt_R.^2./(2.*mu0));...
        rho_R.*vx_R.*vy_R-(Bx.*By_R)./mu0;rho_R.*vx_R.*vz_R-(Bx.*Bz_R)./mu0;...
        vx_R.*By_R-Bx.*vy_R;vx_R.*Bz_R-Bx.*vz_R;
        (E_R+p_R+(Bt_R.^2./(2.*mu0))).*vx_R-((Bt_R.*vt_R)./mu0).*Bx];

    % Add the matrix dissipation term to complete the Roe flux.
    Flux_values = (F_L + F_R  - R*(ev.*CV))/2;
% end
end

% ROE flux
function Flux_values = HLLflux(Q_L,Q_R,gamma)
    % Permeability of free space
    mu0 = 1;
    %mu0 = 4*pi*1e-07;
    Bx = 0.75;
    % Left state
    rho_L = Q_L(1);
    vx_L = Q_L(2)./rho_L;
    vy_L = Q_L(3)./rho_L;
    vz_L = Q_L(4)./rho_L;
    %Bx_L = Q_L(5)./rho_L;
    By_L = Q_L(5);
    Bz_L = Q_L(6);
    E_L = Q_L(7);
    % Total velocity
    vt_L = sqrt((vx_L).^2+(vy_L).^2+(vz_L).^2);
    % Total Magnetic field
    Bt_L = sqrt((Bx).^2+(By_L).^2+(Bz_L).^2);
    B_perp_L = sqrt((By_L).^2+(Bz_L).^2);
    %p_L = (gamma-1)*(Q_L(7) - rho_L*vt_L*vt_L/2 - Bt_L/(2*mu0));
    %p_L = 2/3.*(E_L-0.5.*rho_L.*vt_L.*vt_L-Bt_L.^2/(2.*mu0));
    p_L= 2/3.*(E_L-0.5.*rho_L.*vt_L.*vt_L-Bt_L.^2/(2.*mu0));
    %p_L = rho_L.*E_L+(Bt_L.^2/(2.*mu0));
    a_L = sqrt((gamma*p_L)./rho_L);
    H_L = ( Q_L(7) + p_L );    
    
    % Right state
    rho_R = Q_R(1);
    vx_R = Q_R(2)./rho_R;
    vy_R = Q_R(3)./rho_R;
    vz_R = Q_R(4)./rho_R;
    %Bx_R = Q_R(5)./rho_R;
    By_R = Q_R(5);
    Bz_R = Q_R(6);
    E_R = Q_R(7);
    % Total velocity
    vt_R = sqrt((vx_R).^2+(vy_R).^2+(vz_R).^2);
    % Total Magnetic field
    Bt_R = sqrt((Bx).^2+(By_R).^2+(Bz_R).^2);
    B_perp_R = sqrt((By_R).^2+(Bz_R).^2);
    %p_R = (gamma-1)*(Q_R(7) - rho_R*vt_R*vt_R/2 - Bt_R/(2*mu0));
    p_R = 2/3.*(E_R-0.5.*rho_R.*vt_R.*vt_R-Bt_R.^2/(2.*mu0));
    %p_R = 2/3.*(E_R-0.5.*rho_R.*vt_R.*vt_R-Bt_R.^2/(2.*mu0));
    %p_R = rho_R.*E_R+(Bt_R.^2/(2.*mu0));
    a_R = sqrt((gamma*p_R)./rho_R);
    H_R = ( Q_R(7) + p_R );
    
    % First compute the Roe Averages
    RT = sqrt(rho_R/rho_L);
    rho = RT*rho_L;
    u = (vt_L+RT*vt_R)/(1+RT);
    H = (H_L+RT*H_R)/(1+RT);
    a = sqrt((gamma-1)*(H-u*u/2) );
    
    B_perp_L = sqrt(Bx.^2+By_L.^2);
    B_perp_R = sqrt(Bx.^2+By_R.^2);
    
    Betay_L = By_L/B_perp_L;
    Betaz_L = Bz_L/B_perp_R;
    
    % Alfven speed
    ca_L = sqrt(Bt_L.^2/(4*pi*rho_L));
    ca_R = sqrt(Bt_R.^2/(4*pi*rho_R));
    cax_L = sqrt(Bx.^2/(4*pi*rho_L));
    cax_R = sqrt(Bx.^2/(4*pi*rho_R));
    
    % fast magnetosonic speed
    cf_L = sqrt(1/2.*(a_L.^2+ca_L.^2)+1/2.*sqrt((a_L.^2+ca_L.^2).^2-4*a_L.^2*cax_L.^2));
    cf_R = sqrt(1/2.*(a_R.^2+ca_R.^2)+1/2.*sqrt((a_R.^2+ca_R.^2).^2-4*a_R.^2*cax_R.^2));
    
    % slow magnetosonic speed
    cs_L = sqrt(1/2.*(a_L.^2+ca_L.^2)-1/2.*sqrt((a_L.^2+ca_L.^2).^2-4*a_L.^2*cax_L.^2));
    cs_R = sqrt(1/2.*(a_R.^2+ca_R.^2)-1/2.*sqrt((a_R.^2+ca_R.^2).^2-4*a_R.^2*cax_R.^2));
    
    
    cf_T = sqrt(cf_R/cf_L);
    cf = cf_T*cf_L;
    cs_T = sqrt(cs_R/cs_L);
    cs = cs_T*cs_L;
    ca_T = sqrt(ca_R/ca_L);
    ca = ca_T*ca_L;   
    
     ev = [ abs(u-cf);abs(u-ca);abs(u-cs);abs(u);abs(u+cs);abs(u+ca);abs(u+cf)];
     
    SL = min(vt_L,vt_R)-max(cf_L,cf_R);

    SR = max(vt_L,vt_R)+max(cf_L,cf_R);
    
   
    F_L=[rho_L.*vt_L;rho_L.*vt_L.^2+(Bx.^2./mu0)+p_L+(Bt_L.^2/(2.*mu0));...
         rho_L.*vx_L.*vy_L-(Bx.*By_L)./mu0;rho_L.*vx_L.*vz_L-(Bx.*Bz_L)./mu0;...
         vx_L.*By_L-Bx.*vy_L;vx_L.*Bz_L-Bx.*vz_L;...
         (E_L+p_L+(Bt_L.^2./(2.*mu0))).*vt_L-((Bt_L.*vt_L)./mu0).*Bx];
        
    F_R=[rho_R.*vt_R;rho_R.*vt_R.^2+(Bx.^2./mu0)+p_R+(Bt_R.^2./(2.*mu0));...
        rho_R.*vx_R.*vy_R-(Bx.*By_R)./mu0;rho_R.*vx_R.*vz_R-(Bx.*Bz_R)./mu0;...
        vx_R.*By_R-Bx.*vy_R;vx_R.*Bz_R-Bx.*vz_R;
        (E_R+p_R+(Bt_R.^2./(2.*mu0))).*vt_R-((Bt_R.*vt_R)./mu0).*Bx];
    
    
    if SL<0
        Flux_values = F_L;
    elseif SR>0
        Flux_values = F_R;
    elseif SL<0 && SR>0
        Flux_values = (SR.*F_L-SL.*F_R+SL.*SR.*(Q_R-Q_L))./(SR-SL);
    end
%     % Add the matrix dissipation term to complete the Roe flux.
%     Flux_values = (F_L + F_R)/2  - max(ev).*((Q_R-Q_L)/2);

end



% HLLE Flux
function HLLE = HLLEflux(Q_L,Q_R,gamma)
    % Compute HLLE flux
    % Permeability of free space
    mu0 = 4*pi*1e-07;
    % Left state
    rho_L = Q_L(1);
%     if isnan(rho_L)
%         rho_L = 1e-07;
%     end
%     if rho_L==0
%         rho_L = 1e-07;
%     end
    
    if rho_L~=0 && ~isnan(rho_L)
    vx_L = Q_L(2)./rho_L;
    vy_L = Q_L(3)./rho_L;
    vz_L = Q_L(4)./rho_L;
    Bx_L = 0.75;
    By_L = Q_L(5)./rho_L;
    Bz_L = Q_L(6)./rho_L;
    E_L = Q_L(7)./rho_L;
    % Total velocity
    vt_L = sqrt((vx_L).^2+(vy_L).^2+(vz_L).^2);
    % Total Magnetic field
    Bt_L = sqrt((Bx_L).^2+(By_L).^2+(Bz_L).^2);
    %p_L = (gamma-1)*(Q_L(7) - rho_L*vt_L*vt_L/2 - Bt_L/(2*mu0));
    p_L = 2/3.*(E_L-0.5.*rho_L.*vt_L.*vt_L-Bt_L.^2/(2.*mu0));
    a_L = sqrt(gamma*p_L/rho_L);
    end
    % Right state
    rho_R = Q_R(1);
%     if isnan(rho_R)
%         rho_R = 1e-07;
%     end
%     if rho_R==0
%         rho_R = 1e-07;
%     end
    if rho_R~=0 && ~isnan(rho_L)
    vx_R = Q_R(2)./rho_R;
    vy_R = Q_R(3)./rho_R;
    vz_R = Q_R(4)./rho_R;
    Bx_R = 0.75;
    By_R = Q_R(5)./rho_R;
    Bz_R = Q_R(6)./rho_R;
    E_R = Q_R(7)./rho_R;
    % Total velocity
    vt_R = sqrt((vx_R).^2+(vy_R).^2+(vz_R).^2);
    % Total Magnetic field
    Bt_R = sqrt((Bx_R).^2+(By_R).^2+(Bz_R).^2);
    %p_R = (gamma-1)*(Q_R(7) - rho_R*vt_R*vt_R/2 - Bt_R/(2*mu0));
    p_R = 2/3.*(E_R-0.5.*rho_R.*vt_R.*vt_R-Bt_R.^2/(2.*mu0));
    a_R = sqrt(gamma*p_R/rho_R);
    end
    
    
%     if isnan(rho_L)
%         rho_L = 1e-07;
%     end
%     if isnan(rho_R)
%         rho_R = 1e-07;
%     end
%     if isnan(vx_L)
%         vx_L = 1e-07;
%     end
%     if isnan(vx_R)
%         vx_R = 1e-07;
%     end
%     if isnan(vy_L)
%         vy_L = 1e-07;
%     end
%     if isnan(vy_R)
%         vy_R = 1e-07;
%     end
%     if isnan(vz_L)
%         vz_L = 1e-07;
%     end
%     if isnan(vz_R)
%         vy_L = 1e-07;
%     end
%     
%     if isnan(Bx_L)
%         Bx_L = 1e-07;
%     end
%     
%     if isnan(By_L)
%         By_L = 1e-07;
%     end
%     
%     if isnan(Bz_L)
%         Bz_L = 1e-07;
%     end
%     
%     if isnan(Bx_R)
%         Bx_R = 1e-07;
%     end
%     
%     if isnan(By_R)
%         By_R = 1e-07;
%     end
%     
%     if isnan(Bz_R)
%         Bz_R = 1e-07;
%     end
%     
%     if isnan(vt_L)
%         vt_L = 1e-07;
%     end
%     if isnan(vt_R)
%         vt_R = 1e-07;
%     end
%     if isnan(Bt_L)
%         Bt_L = 1e-07;
%     end
%     if isnan(Bt_R)
%         Bt_R = 1e-07;
%     end
%     if isnan(E_L)
%         E_L = 1e-07;
%     end
%     if isnan(p_L)
%         p_L = 1e-07;
%     end
%     if isnan(p_R)
%         p_R = 1e-07;
%     end
%     
%     if isnan(E_R)
%         E_R = 1e-07;
%     end
    
%if rho_L~=0 && rho_R ~=0 && a_L~=0 && Bt_L~=0 && Bt_R ~=0
    %Alfven speed
    ca_L = sqrt(Bt_L.^2/(4*pi*rho_L));
    ca_R = sqrt(Bt_R.^2/(4*pi*rho_R));
    cax_L = sqrt(Bx_L.^2/(4*pi*rho_L));
    cax_R = sqrt(Bx_R.^2/(4*pi*rho_R));
%end    
    %if ca_L~=0 && ca_R~=0 && a_R~=0 && a_L~=0
    % fast magnetosonic speed
    cf_L = sqrt(1/2.*(a_L.^2+ca_L.^2)+(1/2).*sqrt((a_L.^2+ca_L.^2).^2-4.*a_L.^2.*cax_L.^2));
    cf_R = sqrt(1/2.*(a_R.^2+ca_R.^2)+(1/2).*sqrt((a_R.^2+ca_R.^2).^2-4.*a_R.^2.*cax_R.^2));
%     disp('new');
%     disp(a_R);
%     disp(a_L);
%     disp(ca_L);
%     disp(ca_R);
%     disp(cax_R);
    
    %end
    
%     if isnan(cf_L)
%         cf_L = 1e-07;
%     end
%     if isnan(cf_R)
%         cf_R = 1e-07;
%     end
%     
%     if isnan(vt_L)
%         vt_L = 1e-07;
%     end
%     if isnan(vt_R)
%         vt_R = 1e-07;
%     end
    
    
    SL = min(vt_L,vt_R)-max(cf_L,cf_R);

    SR = max(vt_L,vt_R)+max(cf_L,cf_R);
    disp(SL);
       
    % Flux values
    FL = [rho_L.*vx_L; rho_L.*(rho_L.*vx_L.^2+(Bx_L.^2/mu0)+p_L+(Bt_L.^2/(2*mu0)));...
        rho_L.*(rho_L.*vx_L.*vy_L-Bx_L.*By_L/mu0);rho_L.*(rho_L.*vx_L.*vz_L-Bx_L.*Bz_L/mu0);...
        vx_L.*By_L-Bx_L.*vy_L;vx_L.*Bz_L-Bx_L.*vz_L;(E_L+p_L+Bt_L.^2/(2*mu0)).*vx_L-((Bt_L.*vt_L)/mu0).*Bx_L];
    
    FR = [rho_R.*vx_R; rho_R.*(rho_R.*vx_R.^2+(Bx_R.^2/mu0)+p_R+(Bt_R.^2/(2*mu0)));...
        rho_R.*(rho_R.*vx_R.*vy_R-Bx_R.*By_R/mu0);rho_R.*(rho_R.*vx_R.*vz_R-Bx_R.*Bz_R/mu0);...
        vx_R.*By_R-Bx_R.*vy_R;vx_R.*Bz_R-Bx_R.*vz_R;(E_R+p_R+Bt_R.^2/(2*mu0)).*vx_R-((Bt_R.*vt_R)/mu0).*Bx_R];
    
    % Compute the HLL flux.
    if SL >= 0  % Right-going supersonic flow
        HLLE = FL;
    elseif (SL <= 0) && (SR >= 0) % Subsonic flow
        HLLE = ( SR*FL - SL*FR + SL*SR*(Q_R-Q_L) )/(SR-SL);

    elseif  SR <= 0 % Left-going supersonic flow
        HLLE = FR;
    end
    
    % Special case Lax-Friedrich flux
    if SR==-SL
        HLLE = (FL+FR)/2-SR*(Q_R-Q_L)/2;
    end
%     if isnan(SL) || isnan(SR)
%         HLLE = 0;
%     end
    
end

