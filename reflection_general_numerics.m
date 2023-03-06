clear;
close all

% % % % case 0 ( symmetry)
% dd1 = 1;      dd2 = 3;
% cc1 = 1;      cc2 = 2;
% mm1 = 1;      mm2 = 4;

% dd1 = 3;      dd2 = 1; % in opposite direction 
% cc1 = 2;      cc2 = 1;
% mm1 = 4;      mm2 = 1;

% % % case 1
% dd1 = 1;      dd2 = 1;
% cc1 = 1;      cc2 = 1;
% mm1 = 1;      mm2 = 2;

% % % % case 2
dd1 = 0.5;    dd2 = 1;
cc1 = 2;      cc2 = 1;
mm1 = 3;      mm2 = 1;
cc1_mod = cc2;

% % test !!!!
% dd1 = 0.5;    dd2 = 1;
% cc1 = 2;      cc2 = 1;
% mm1 = 3;      mm2 = 1;
% test !!!!

% % % % case 3
% dd1 = 2;      dd2 = 0;
% cc1 = 1;      cc2 = 1;
% mm1 = 3;      mm2 = 2;


% % % transparency
% dd1 = 3;
% mm1 = 3;
% dd2 = 0;
% mm2 = 1;
% tild_om_t = 0.25;
% cc1 = 1/tild_om_t * 0.25*(mm2*dd1-dd2*mm1)/(mm1 - mm2);
% cc2 = cc1;

% dd1 = 1;      dd2 = 1;
% cc1 = 1;      cc2 = 1;
% mm1 = 1;      mm2 = 1;


beta_   = 0.02; %0.02-good % ICs - u_nn = exp(-beta^2/2*((nn-nn0)*a-vg*t)^2)*sin(k*nn-Omega*t) 
nn0     = -4/beta_; % width of the the wave packet is ~3-4 * 1/beta
aa        = 1;
Nom       = 20;
write_gif = 0;
N         = 1000; 
ii_step = 1000;

% -------------------------------------------------------------------------
omega_min1 = sqrt((dd1)/mm1);                 omega_min2 = sqrt((dd2)/mm2);
omega_max1 = sqrt((4*cc1 + dd1)/mm1);         omega_max2 = sqrt((4*cc2 + dd2)/mm2);
omega_min  = max(omega_min1, omega_min2);     omega_max  = min(omega_max1, omega_max2);

dt        = 0.05/omega_max;

Omega    = zeros(Nom,1);    tilde_Omega_2 = zeros(Nom,1);
v1_exact = zeros(Nom,1);    v2_exact = zeros(Nom,1);
eps_num = zeros(Nom,1);
T_an    = zeros(Nom,1);     T_Sim   = zeros(Nom,1);     T_num   = zeros(Nom,1);
T_an_new = zeros(Nom,1);

filename = 'reflection_gif.gif';
if write_gif == 2
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
end

for kk = 1:Nom
%     Omega(kk) = sqrt((dd2-dd1)/(mm2-mm1));%omega_min + kk*(omega_max-omega_min)/(Nom + 1);
    Omega(kk) = omega_min + kk*(omega_max-omega_min)/(Nom + 1);

    theta     = acos(mm1/(2*cc1) * ((2*cc1 + dd1)/mm1 - Omega(kk)^2));    
  % ----------   analytics   -----------------
    v1(kk)   = 1/(2*Omega(kk)) *sqrt((Omega(kk)^2-dd1/mm1)*(4*cc1/mm1 + dd1/mm1 - Omega(kk)^2));
    v2(kk)   = 1/(2*Omega(kk)) *sqrt((Omega(kk)^2-dd2/mm2)*(4*cc2/mm2 + dd2/mm2 - Omega(kk)^2));
    
    G_G      =   aa^2/(2*v2(kk)) * ...
                  (...
                  ((cc2 - cc1)/mm2 + (mm1-mm2)*cc2/(mm1*mm2)) *(1 - dd2/(mm2*Omega(kk)^2))...
                  +...
                   0.5/(Omega(kk)^2) *(dd1/mm1 - dd2/mm2)* (2*cc1/mm2 + dd2/mm2 - Omega(kk)^2)...
                  );
              
    % CHECK !! ! !  !!    
%     cc12 = cc1_mod;
%     GG1 = aa^2/(2*v2(kk)*Omega(kk)^2) * ((...
%                             cc2*(2*mm1 - mm2) - mm1*cc12)/(mm1*mm2)*(Omega(kk)^2 - dd2/mm2)...
%                             +...
%                             0.5*((2*cc12 + dd2)/mm2 - Omega(kk)^2)*(dd1/mm1 - dd2/mm2));
%                         
%     HH = 0.5 * aa^2 * (cc1 - cc12) * (Omega(kk)^2 - dd1/mm1) / (mm1 * Omega(kk)^2 * v1(kk));
%     T_an_new(kk) = (2*v1(kk)) / (v1(kk) + v2(kk) - GG1 + HH);

   
    E1_E2    = 1/(2*v1(kk)) * (v2(kk) - v1(kk) - G_G);
%     E1_E2    = 1/(2*v1(kk)) * ( v2(kk) - v1(kk) - aa^2/(2*v2(kk)) * ...
%                   (...
%                   ((cc2 - cc1)/mm2 + (mm1-mm2)*cc2/(mm1*mm2)) *(1 - dd2/(mm2*Omega(kk)^2))...
%                   +...
%                    0.5/(Omega(kk)^2) *(dd1/mm1 - dd2/mm2)* (2*cc1/mm2 + dd2/mm2 - Omega(kk)^2)...
%                   ));
%       
              
              
    T_an(kk) = 1/(E1_E2 + 1);    
    T_Sim(kk) = 4*mm1*mm2*v1(kk)*v2(kk) / (((mm1-mm2)*Omega(kk)^2 + dd2-dd1)^2/(4*Omega(kk)^2) +(mm1*v1(kk) + mm2*v2(kk))^2);
% ----------------------------------------------
    % ------------   initial conditions    ----------------
    vp        = zeros(3*N + 1,1);    up    = zeros(3*N + 1,1);
    F2        =  zeros(3*N + 1,1);
%     WP_length = round(N_waves * 2 * pi/theta);
    
    for jj = 1 : N
        nn  = jj - N - 1;
        k_ = theta;
        up(jj) =  exp(-0.5*beta_^2*(nn-nn0)^2) * sin(k_*nn);
        vp(jj) =  exp(-0.5*beta_^2*(nn-nn0)^2) *(-Omega(kk)* cos(k_*nn)+...
                 + beta_^2*(nn-nn0) * v1(kk) * sin(k_*nn));
%         kx = theta * (jj - N) + pi * N_waves;
%         up(jj) =  exp(-beta*kx^2) * sin(kx);
%         vp(jj) = -Omega(kk) * exp(-beta*kx^2) * (cos(kx) - 2 * beta * kx * sin(kx));
    end
    
    if write_gif == 0 
        figure(100); plot( (N+1:3*N+1) - N-1, up(N+1:3*N+1), '.-', 'Color', [0 0.5 0]);
        hold on
        plot((1:N) - N-1,up(1:N), 'b.-',  [-0.1 +0.1], [-1; 1], '--');
        hold off
        ylim([-1; 1]);
    end
%     pause

    % -----------------------------------------------------      
    tmax = min(20*abs(nn0)/v1(kk), 20*abs(nn0));
    Ndt  = round(tmax/dt);
    t    = ((1:Ndt)*dt)' - dt;
       
    Kin_1 = zeros(Ndt,1);      Kin_2 = zeros(Ndt,1);     
    Pot_1 = zeros(Ndt,1);      Pot_2 = zeros(Ndt,1);       
    % ----------------------------------------------------
    
    diff_E2 = zeros(Ndt, 1);
    Lagrangian = zeros(Ndt,1);
    
    ii = 1;
    stop_flag = 0;
    while stop_flag == 0 && ii < Ndt       
        vp(1) = vp(1) + (cc1/mm1) * (up(2) -  up(1)) * dt;
        for jj=2:N-1
            vp(jj) = vp(jj) - (dd1/mm1)*up(jj)*dt + (cc1/mm1) *(up(jj+1) - 2 * up(jj) + up(jj-1)) * dt;
        end     

        % в аналитике это нулевая частица соответствует частице N+1 в моделировании:
        vp(N)   = vp(N)   - (dd1/mm1)*up(N)*dt   + (cc1_mod*(up(N+1) - up(N)) + cc1*(up(N-1)-up(N)))     * dt/mm1;
        vp(N+1) = vp(N+1) - (dd2/mm2)*up(N+1)*dt + (cc2*(up(N+2) -  up(N+1))  + cc1_mod*(up(N)- up(N+1))) * dt/mm2;
        
        for jj=N+2:3*N
            vp(jj) = vp(jj) - (dd2/mm2) * up(jj) * dt + (cc2/mm2) * (up(jj+1) - 2 * up(jj) + up(jj-1)) * dt;        end
        for jj=1:3*N
            up(jj) = up(jj) + vp(jj) * dt;
        end
        
% ------------------------    Energies    --------------------------------------
        for jj=1:N-1
            Pot_1(ii) = Pot_1(ii) + 0.5*dd1*up(jj)^2 + 0.5 * cc1 * (up(jj+1) - up(jj))^2;        
        end
        Pot_1(N) = Pot_1(N) + 0.5*dd1*up(N)^2 + 0.5 * cc1_mod * (up(N+1) - up(N))^2;
        for jj=N+1:3*N-1
            Pot_2(ii) = Pot_2(ii) + 0.5*dd2*up(jj)^2 + 0.5 * cc2 * (up(jj+1) - up(jj))^2;
        end
        Kin_1(ii)  = 0.5 * mm1 * sum(vp(1:N).^2);
        Kin_2(ii)  = 0.5 * mm2 * sum(vp(N+1:3*N).^2);
        
        if (mod(ii,ii_step) == 0)  disp(ii);
        end
                    
        if (mod(ii,ii_step) == 0 && write_gif == 1)
%             width_ = 4*WP_length;
            plot( (N+1:3*N+1) - N, vp(N+1:3*N+1)/Omega(kk), '.-', 'Color', [0 0.5 0]);
            hold on
            plot((1:N) - N, vp(1:N)/Omega(kk), 'b.-',  [-0.1 +0.1], [-1; 1], '--');
            hold off
            xlim([ -N; 2*N]);
            ylim([-1; 1]);
            ylabel('v_{n}/(a \Omega)')
            title(['\Omega/\Omega_{max} = ', num2str(Omega(kk)/omega_max), '   T = ', num2str(T_an(kk)),...
                '   d_1=', num2str(dd1), '  c_1=',num2str(cc1), '  m_1=', num2str(mm1),...
                '  d_2=', num2str(dd2), '  c_2=', num2str(cc2), '  m_2=', num2str(mm2) ]);
            set(gcf, 'color', [1 1 1]);
            daspect([900 1 1]);
            drawnow
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ii == ii_step && kk == 1
              imwrite(imind,cm,filename,'gif','DelayTime',0.001, 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','DelayTime',0.001, 'WriteMode','append');
            end
        end        
        % -------------------------------------------------------------
                
        if ii > 4*abs(nn0)/(v1(kk)*dt)
            stop_flag = 1;
        end
        ii = ii + 1;
    end

   eps_num(kk) = (Kin_1(ii-1) + Pot_1(ii-1))/(Kin_2(ii-1) + Pot_2(ii-1));
   T_num(kk)   = 1/(1 + eps_num(kk)); 
   
   tilde_Omega_2(kk) = (Omega(kk)^2 - omega_min^2)/(omega_max^2 - omega_min^2); 

   figure(5);  plot(t(1:ii), (Kin_1(1:ii) + Pot_1(1:ii))/(Kin_1(1) + Pot_1(1)), t(1:ii), (Kin_2(1:ii) + Pot_2(1:ii))/(Kin_1(1) + Pot_1(1)));
   title(['T = ', num2str(T_num(kk))]);
   ylim([0; 1]);
   set(gcf, 'color', [1 1 1]);
   figure(55); plot(up);
%    figure(56); plot(1:3*N+1, F2, N, F2(N),'ro', N+1, F2(N+1), 'g+'); 
%    figure(57); plot(Lagrangian); title('Lagrangian of particle 0');
%    pause;
end

figure(6);   plot(tilde_Omega_2, T_num, 'r.', tilde_Omega_2, T_an, 'm.', tilde_Omega_2, T_Sim, 'b-');%, Lambda/om_max_2, eps_cont*ones(Nom, 1));
xlabel('\Omega/\omega_{max}');
ylabel('T = E_2/E');
title(['d_1=', num2str(dd1), '  c_1=',num2str(cc1), '  m_1=', num2str(mm1),...
    '  d_2=', num2str(dd2), '  c_2=', num2str(cc2), '  m_2=', num2str(mm2),... 
    '  \max T_{num} = ', num2str(max(T_num)), '  \max T_{ab} = ', num2str(max(T_an))]);
ylim([0 1]);


for ii=1:length(T_an)
    rel_er(ii) = (T_an(ii) - T_num(ii))/T_an(ii);
end


% figure(7);   plot(tilde_Omega_2, rel_er, 'r.');%, Lambda/om_max_2, eps_cont*ones(Nom, 1));
% xlabel('\Omega/\omega_{max}');
% ylabel('T = E_2/E');
% title(['d_1=', num2str(dd1), '  c_1=',num2str(cc1), '  m_1=', num2str(mm1),...
%     '  d_2=', num2str(dd2), '  c_2=', num2str(cc2), '  m_2=', num2str(mm2),... 
%     '  \max T_{num} = ', num2str(max(T_num)), '  \max T_{ab} = ', num2str(max(T_an))]);

% 
% 
% 
% figure(7);   plot(Omega/omega_max, T_num, 'r.', Omega/omega_max, T_an, 'b.', Omega/omega_max, T_inclusion, 'mo');%, Lambda/om_max_2, eps_cont*ones(Nom, 1));
% xlabel('\Omega/\omega_{max}');
% ylabel('T = E_2/E');
% title(['d_1=', num2str(dd1), '  c_1=',num2str(cc1), '  m_1=', num2str(mm1),...
%     '  d_2=', num2str(dd2), '  c_2=', num2str(cc2), '  m_2=', num2str(mm2)]);
% ylim([0 1]);
% % figure(7);   plot(Omega/omega_max, 100*(T_an - T_num).*T_num.^(-1), 'r.', Omega/omega_max, 100*(T_Sim - T_num).*T_num.^(-1), 'g.'); 
% % title('Error, %');
% 
% 
% % figure(8); plot(Omega/omega_max, mm1*v1, Omega/omega_max, mm2*v2);