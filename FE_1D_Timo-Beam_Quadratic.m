%clear all; close all;

% Setting parameters for solid
E = 200e9;          % Youngs Mod
t=0.1;              % Thickness of beam
w=0.2;              % Width of beam
A = t*w;            % Area of beam
k = 5/6;            % Timoshenko shear coefficient
Gg = 79.3e9;        % Shear modulus of steel
I = t^3*w/12;       % Moment of inertia
Le = 1;             % Length of bar
P=0;               % Force on end of beam
rho = -1;

theta_bc = 0;
u_bc = 0;
mbc = 0;            %/(E*I)?
sbc = P;

% Setting parameters for FEM
nel = 2;               % number elements
dof = nel*2+1;         % number of nodes/matrix size
x = linspace(0,Le,dof); % nodal x-points
MCA = [1:dof-2; 2:dof-1; 3:dof]';
y0 = zeros(1,3);
el_l = Le/nel;
eps_inv = 1/.000000001;       % k*Gg*A*el_l^2/(E*I)

hh = waitbar(0,'Building Stiffness and Mass Matrices. Please Wait...');

GP = [-1/sqrt(3)  1/sqrt(3)]; % 2 point***********
W = [1  1];    % 2 point ************
% GP = [-.7745966692 0 .7745966692]; % 3 point
% W = [.5555555556 .8888888889 .5555555556];    % 3 point
% GP = [-.8611363116 -.3399810436 .8611363116 .3399810436]; % 4 point
% W = [.3478548451 .6521451549 .3478548451 .6521451549];    % 4 point

K = zeros(dof);
L = zeros(dof);
M = zeros(dof);

F = zeros(dof,1);
G = zeros(dof,1);
elx = zeros(nel,3);

set(gcf, 'Visible', 'off')                  % make plotting invisible
figure(1)
for el = 1:nel % Loop through each element
    
    elx(el,1) = x(1,el*2-1);        % nodal coord array x1 point
    elx(el,2) = x(1,el*2);      % nodal coord array x2 point
    elx(el,3) = x(1,el*2+1);      % nodal coord array x2 point
    
    plot (elx(el,1:3), y0,'Color','k') % plot undeformed mesh
    hold on                                     
    
    %el_l = elx(el,2)-elx(el,1);                             % element length

    K_e = zeros(3);
    L_e = zeros(3);
    M_e = zeros(3);
    
    F_e = zeros(3,1);
    G_e = zeros(3,1);
    for i=1:max(size(GP))

            psi = GP(i);
            
             % Parent element shape func matrix
            N = [.5*(psi-1)*psi -(psi-1)*(1+psi) .5*(psi+1)*psi];
            detJ = el_l/2;         % Jacobian inverse
            invdetJ = detJ^-1;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Parent element grad shape func matrix
            B = [psi-.5 -2*psi .5+psi]*invdetJ;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            K_e = K_e + W(i)*(B'*B)*detJ;
            L_e = L_e + W(i)*(N'*B)*detJ;
            M_e = M_e + W(i)*(N'*N)*detJ;
            F_e = F_e + W(i)*N'*rho*detJ;
            
    end
        % Implementing Neumann BCs into force vectors
        if el == nel
            psi = 1;
            % Parent element shape func matrix
            N = [.5*(psi-1)*psi -(psi-1)*(1+psi) .5*(psi+1)*psi];
            F_e = F_e + N'*sbc;
            G_e = G_e + N'*mbc;
        end
    
    K(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) = K(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) + K_e;
    L(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) = L(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) + L_e;
    M(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) = M(MCA(el*2-1,1:3),MCA(el*2-1,1:3)) + M_e;
    
    F(MCA(el*2-1,1:3),1) = F(MCA(el*2-1,1:3),1) + F_e;
    G(MCA(el*2-1,1:3),1) = G(MCA(el*2-1,1:3),1) + G_e;
    waitbar(el/nel);
end
K11 = K + eps_inv*M;
K12 = -eps_inv*L;
K22 =  eps_inv*K;
K_f = [K11 , K12 ; K12' , K22];

F(1) = 0;
F(dof) = 0;

f_f = [G ; F];

K_f(1,:) = 0;
K_f(:,1) = 0;
K_f(1,1) = 1;

K_f(dof+1,:) = 0;
K_f(:,dof+1) = 0;
K_f(dof+1,dof+1) = 1;

K_f(dof,:) = 0;
K_f(:,dof) = 0;
K_f(dof,dof) = 1;

K_f(dof+dof,:) = 0;
K_f(:,dof+dof) = 0;
K_f(dof+dof,dof+dof) = 1;
% K_f(1,dof+1) = 1;
% K_f(dof+1,1) = 1;

f_f = f_f - K_f(:,1)*theta_bc - K_f(:,dof+1)*u_bc;

d = K_f\f_f;
format long
min(d(dof:dof+dof))
close(hh);

hold on
for el = 1:nel
    plot (elx(el,1:3), d(dof+el*2-1:dof+el*2+1),'Color','r') % plot undeformed mesh
    hold on    
end
set(gcf, 'Visible', 'on')
% Create xlabel
xlabel({'x (m)'},'FontWeight','demi','FontSize',14);
% Create ylabel
ylabel({'y (m)'},'FontWeight','demi'...
    ,'FontSize',14);
print('exact', '-dpng', '-r600');
%legend('\epsilon = 0.1', '\epsilon = 0.01', '\epsilon = 0.001')