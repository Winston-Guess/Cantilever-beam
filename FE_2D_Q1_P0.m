%clear all; close all;

% Setting parameters for solid
E = 1500;          % Youngs Mod
v = 0.499;
height = 2;
Length = 10;              % length of bar
% 
lamb = E*v/((1+v)*(1-2*v));
mu = E/(2*(1+v));
Dmu = [2*mu 0 0;0 2*mu 0;0 0 mu];

f_ = 10;            % Distributed load max (f bar)
%disp(D);

nel_y = 2;   nnp_y = nel_y + 1;     % number of elem and nodes in y-dir
nel_x = 4;  nnp_x = nel_x + 1;     % number of elem and nodes in x-dir
nel = nel_y * nel_x;                % number total elem
nnp = nnp_y * nnp_x;                % number total nodes
Le = min(height/nel_y,Length/nel_x);         % element length in x-dir
dof = 2*(nel_y+1)*(nel_x+1);      % size of M and K matrices 

NPxy = zeros(2,nnp);            % Nodal point x and y values
ICA = zeros(nel,4);             % Nodal Interconnectivity Matrix
MCA = zeros(nel,8);             % Nodal xy Interconnectivity Matrix
Ess_BC_y = zeros(1,nnp);            % y Essential boundary vetor
Ess_BC_x = zeros(1,nnp);            % x Essential boundary vetor
Nat = zeros(1,nnp);             % Natural boundary vetor

EBC_yCount = 0;
EBC_xCount = 0;

% Algorithm for flagging Nodal Boundary Conditions
for np = 1:nnp
    NPxy(1,np) = EBC_xCount*Length/nel_x;
    if NPxy(1,np) == 10
        Nat(1,np) = 1;
    else if NPxy(1,np) == 0
            Ess_BC_x(1,np) = 1;
        end
    end 
    NPxy(2,np) = (nnp_y-1-EBC_yCount)*height/(nnp_y-1);
    if NPxy(2,np) == 0 && NPxy(1,np) == 0
        Ess_BC_y(1,np) = 1;
    end 
    if EBC_yCount == nnp_y-1
        EBC_yCount = 0;
        EBC_xCount = EBC_xCount+1; 
    else
        EBC_yCount = EBC_yCount+1;
    end
end
colmCor = 0;  % Corrector for connectivity matrices based on column counter
colmCount = 1;   % Counter for the number of columns of elements (el x-dir)
% Algorithm for generating connectivity matrices
for np = 1:nel
    ICA(np,1) = np + (colmCor);
    MCA(np,1) = (np + (colmCor))*2-1;
    MCA(np,2) = (np + (colmCor))*2;
    
    ICA(np,2) = np+1 + (colmCor);
    MCA(np,3) = (np+1 + (colmCor))*2-1;
    MCA(np,4) = (np+1 + (colmCor))*2;
    
    ICA(np,3) = np+1+nnp_y + (colmCor);
    MCA(np,5) = (np+1+nnp_y + (colmCor))*2-1;
    MCA(np,6) = (np+1+nnp_y + (colmCor))*2;
    
    ICA(np,4) = np +nnp_y + (colmCor);
    MCA(np,7) = (np+nnp_y + (colmCor))*2-1;
    MCA(np,8) = (np+nnp_y + (colmCor))*2;
    
    if colmCount == nnp_y-1
        colmCount = 1;
        colmCor = colmCor+1;
    else
        colmCount = colmCount+1;
    end
end
   
elx = zeros(nel,5);
ely = zeros(nel,5);

hh = waitbar(0,'Building Stiffness and Mass Matrices. Please Wait...');

GP = [-1/sqrt(3)  1/sqrt(3)]; % 2 point***********
W = [1  1];    % 2 point ************
% GP = [-.7745966692 0 .7745966692]; % 3 point
% W = [.5555555556 .8888888889 .5555555556];    % 3 point
% GP = [-.8611363116 -.3399810436 .8611363116 .3399810436]; % 4 point
% W = [.3478548451 .6521451549 .3478548451 .6521451549];    % 4 point

    F = zeros(dof,1);
    Kmu = zeros(dof);
    M = eye(nel);
    G = zeros(dof,nel);
set(gcf, 'Visible', 'off')
figure(1)
for el = 1:nel
    for i=1:4
        elx(el,i) = NPxy(1,ICA(el,i));
        ely(el,i) = NPxy(2,ICA(el,i));

        if i==1
            elx(el,5) = NPxy(1,ICA(el,i));
            ely(el,5) = NPxy(2,ICA(el,i));
        end
    end
    set(gcf, 'Visible', 'off')
    plot (elx(el,1:5), ely(el,1:5),'Color','k')
    hold on
    
    h = ely(el,4)-ely(el,3);

    xy_e = transpose(NPxy(1:2,ICA(el,1:4)));

    K_e_mu = zeros(8);
    F_e = zeros(8,1);
    L_e = zeros(1,8);
    
    for i=1:max(size(GP))
        for j=1:max(size(GP))
            % Guass points in this loop
            psi = GP(i);
            eta = GP(j);
            % GN matrix derived form the gradient operator acting on the
            % shape function matrix
            GN = .25*[eta-1 1-eta 1+eta -eta-1
                      psi-1 -psi-1 1+psi 1-psi];
            J = GN*xy_e;            % Jacobian
            detJ = det(J);          % Determinant jacobian        
         
            BB = J\GN;% compute the derivative of the shape functions

            B1x = BB(1,1); B2x = BB(1,2); B3x = BB(1,3); B4x = BB(1,4);

            B1y = BB(2,1); B2y = BB(2,2); B3y = BB(2,3); B4y = BB(2,4);
        
            B = [ B1x      0     B2x     0      B3x    0      B4x     0  ;
                0     B1y     0     B2y      0     B3y     0      B4y; 
              B1y     B1x    B2y    B2x     B3y    B3x    B4y     B4x];
            
            % Summing up the element stiffness and mass matrices
            K_e_mu = K_e_mu + W(i)*W(j)*(B'*Dmu*B)*detJ;
            
            L_e = L_e + [B1x B1y, B2x B2y, B3x B3y, B4x B4y];

        end
        
        if Nat(ICA(el,4)) == 1  % Checks if element has a natural BC
            eta=1;    % Set eta to one for guass integration along that edge
            % Adjusting shape funtions for building up the element force matrix
            N1 = .25*(1-psi)*(1-eta);
            N2 = .25*(1+psi)*(1-eta);
            N3 = .25*(1+psi)*(1+eta);
            N4 = .25*(1-psi)*(1+eta);

            yN = [ely(el,1) 0 ely(el,2) 0 ely(el,3) 0 ely(el,4) 0];

            N = [N1 0 N2 0 N3 0 N4 0;0 N1 0 N2 0 N3 0 N4];
            
            fyN = -f_*(yN*N'*[1; 0]-1);
            disp(fyN)

            F_e = F_e + W(i)*fyN*N'*[1; 0]*h/2;
            
        end
    end 
    %disp(F_e)
    % Scatter and add the element stiff matrix into the global stiff matrix
    Kmu(MCA(el,1:8),MCA(el,1:8)) = Kmu(MCA(el,1:8),MCA(el,1:8)) + K_e_mu;
    % Scatter and add the element force vector into the global force vector
    G(MCA(el,1:8),el) = G(MCA(el,1:8),el) - L_e';
    % Scatter and add the element forc e vector into the global force vector
    F(MCA(el,1:8),1) = F(MCA(el,1:8),1) + F_e;
    waitbar(el/nel);
end
close(hh);
%disp(F)
% Algorithm: Applies partioning to K matrix at Ess BCs x and y

Kf = Kmu -lamb*G*(M\G');

for np = 1:nnp
    % Checking for true flags in x-dir essential boundary flag vec
    if Ess_BC_x(np) == 1
        Kf(np*2-1,:) = 0;
        Kf(:,np*2-1) = 0;
        Kf(np*2-1,np*2-1) = 1;
    end
    % Checking for true flags in y-dir essential boundary flag vec
    if Ess_BC_y(np) == 1
        Kf(np*2,:) = 0;
        Kf(:,np*2) = 0;
        Kf(np*2,np*2) = 1;
    end
end

%K = sparse(K);
%  F(dof-5) = -10;
%  F(dof-3) = -20;
%  F(dof-1) = -10;
S=1;
Ff = F;
df = Kf\Ff;
df=S*df;
set(gcf, 'Visible', 'off')
for el = 1:nel
    for i=1:4
        elx(el,i) = elx(el,i)+ df(ICA(el,i)*2-1);
        ely(el,i) = ely(el,i) +df(ICA(el,i)*2);

        if i==1
            elx(el,5) = elx(el,5)+df(ICA(el,i)*2-1);
            ely(el,5) = ely(el,5)+df(ICA(el,i)*2);
        end
    end
    set(gcf, 'Visible', 'off')
    plot (elx(el,1:5), ely(el,1:5),'Color','r')
    hold on
end
set(gcf, 'Visible', 'on')
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',14);
% Create ylabel
ylabel({'y (m)'},'FontWeight','normal'...
    ,'FontSize',14);
x= 10;
y=2;
u_xy = 2*f_*(1-v^2)/(E*height)*x*(height/2-y);
disp('Analytical x displacement = ')
disp(u_xy)
disp('FE x displacement = ')
disp(df(13*2-1))
v_xy = f_*(1-v^2)/(E*height)*(x^2+v/(1-v)*y*(y-height));
disp('Analytical y displacement = ')
disp(v_xy)
disp('FE y displacement = ')
disp(df(13*2))