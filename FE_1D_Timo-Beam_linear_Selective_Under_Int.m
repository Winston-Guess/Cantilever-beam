clear all; close all;

% Setting parameters for solid
% E = 200e9;          % Youngs Mod
% t=0.1;              % Thickness of beam
% w=0.2;              % Width of beam
% A = t*w;            % Area of beam
% k = 5/6;            % Timoshenko shear coefficient
% Gg = 79.3e9;        % Shear modulus of steel
% I = t^3*w/12;       % Moment of inertia
Le = 1;               % Length of bar
P=-1;                 % Force on end of beam
rho = -1;             % Uniform dist load

theta_bc = 0;
u_bc = 0;
mbc = 0;            %/(E*I)?
sbc = 0;

% Setting parameters for FEM
nel = 100;               % number elements
dof = nel+1;         % number of nodes/matrix size
x = linspace(0,Le,dof); % nodal x-points
MCA = [1:dof-1; 2:dof]';
y0 = zeros(1,2);
el_l = Le/nel;
eps_inv =  1/.1;     % k*Gg*A*Le^2/(E*I); 

hh = waitbar(0,'Building Stiffness and Mass Matrices. Please Wait...');

GP1 = 0;
W1 = 2;
GP = [-1/sqrt(3)  1/sqrt(3)]; % 2 point***********
W = [1  1];    % 2 point ************
% GP = [-.7745966692 0 .7745966692]; % 3 point
% W = [.5555555556 .8888888889 .5555555556];    % 3 point
% GP = [-.8611363116 -.3399810436 .8611363116 .3399810436]; % 4 point
% W = [.3478548451 .6521451549 .3478548451 .6521451549];    % 4 point

K = zeros(dof);
K_UI = zeros(dof);
L = zeros(dof);
M = zeros(dof);

F = zeros(dof,1);
G = zeros(dof,1);
elx = zeros(nel,2);

% set(gcf, 'Visible', 'off')                  % make plotting invisible
% figure(1)
for el = 1:nel % Loop through each element
    for i=1:2       
        elx(el,1) = x(1,el);        % nodal coord array x1 point
        elx(el,2) = x(1,el+1);      % nodal coord array x2 point
    end
%     plot (elx(el,1:2), y0,'Color','k') % plot undeformed mesh
%     hold on                                     
    
    %el_l = elx(el,2)-elx(el,1);          % element length

    K_e = zeros(2);   
    F_e = zeros(2,1);
    G_e = zeros(2,1);
    for i=1:max(size(GP))

            psi = GP(i);

            N = .5*[ 1-psi 1+psi ]; % Parent element shape func matrix
            detJ = el_l/2;         % Jacobian inverse
            invdetJ = detJ^-1;
            % Parent element grad shape func matrix
            B = .5*[ -1 1 ]*invdetJ; %          
            
            K_e = K_e + W(i)*(B'*B)*detJ;
            F_e = F_e + W(i)*N'*rho*detJ;
            
    end
    
    psi = GP1;
    N = .5*[ 1-psi 1+psi ]; % Parent element shape func matrix
    detJ = el_l/2;         % Jacobian inverse
    invdetJ = detJ^-1;
    % Parent element grad shape func matrix
    B = .5*[ -1 1 ]*invdetJ; %  
    
    K_e_UI = W1*(B'*B)*detJ;
    L_e = W1*(N'*B)*detJ;
    M_e = W1*(N'*N)*detJ;
    
    if el == nel
        psi = -1;
        N = .5*[ 1-psi 1+psi ]; % Parent element shape func matrix
        F_e = F_e + N'*sbc;
        G_e = G_e + N'*mbc;
    end
    
    K(MCA(el,1:2),MCA(el,1:2)) = K(MCA(el,1:2),MCA(el,1:2)) + K_e;
    L(MCA(el,1:2),MCA(el,1:2)) = L(MCA(el,1:2),MCA(el,1:2)) + L_e;
    M(MCA(el,1:2),MCA(el,1:2)) = M(MCA(el,1:2),MCA(el,1:2)) + M_e;
    K_UI(MCA(el,1:2),MCA(el,1:2)) = K_UI(MCA(el,1:2),MCA(el,1:2)) + K_e_UI;
    
    F(MCA(el,1:2),1) = F(MCA(el,1:2),1) + F_e;
    G(MCA(el,1:2),1) = G(MCA(el,1:2),1) + G_e;
    waitbar(el/nel);
end
K11 = K + eps_inv*M;
K12 = -eps_inv*L;
K22 =  eps_inv*K_UI;
K_f = [K11 , K12 ; K12' , K22];

f_f = [G ; F];

f_f = f_f - K_f(:,1)*theta_bc - K_f(:,dof+1)*u_bc-...
        K_f(:,dof)*theta_bc - K_f(:,dof+dof)*u_bc;

%BC at node 1
K_f(1,:) = 0;
K_f(dof+1,:) = 0;
K_f(:,1) = 0;
K_f(:,dof+1) = 0;
K_f(1,1) = 1;
K_f(dof+1,dof+1) = 1;
%BC at node n
K_f(dof,:) = 0;
K_f(dof+dof,:) = 0;
K_f(:,dof) = 0;
K_f(:,dof+dof) = 0;
K_f(dof,dof) = 1;
K_f(dof+dof,dof+dof) = 1;
% K_f(1,dof+1) = 1;
% K_f(dof+1,1) = 1;
f_f(dof+1) = 0;
f_f(dof+dof) = 0;

d = K_f\f_f;
format long
min(d(dof:dof+dof))
close(hh);

figure(1)
for el = 1:nel
    plot (elx(el,1:2), d(dof+el:dof+el+1),'Color','k')% plot undef mesh
    hold on    
end
% set(gcf, 'Visible', 'on')
figure(2)
for el = 1:nel
    plot (elx(el,1:2), d(el:el+1),'Color','k')% plot undef mesh
    hold on    
end

