clear all ; close all ; clc ;

% Set the number of elements
n_el = 16; n_np = 15;
dof = 2 * n_np;
f_ = 1;
height = 2;
Length = 10;              % length of bar
% Set the properties of the material and problem
E = 1500;          % Youngs Mod
v = .4;
% 
lamb = E*v/((1+v)*(1-2*v));
mu = E/(2*(1+v));
Dmu = [2*mu 0 0;0 2*mu 0;0 0 mu];

ICA = [1 5 4            % Interconnectivity array
       1 2 5
       2 6 5
       2 3 6
       4 8 7
       4 5 8
       5 9 8
       5 6 9
       7 11 10
       7 8 11
       8 12 11
       8 9 12
       10 14 13
       10 11 14
       11 15 14
       11 12 15];

dofICA = [ICA*2-1,ICA*2];
   
GA = [ 0  0  0  2.5  2.5  2.5  5  5  5  7.5  7.5  7.5  10  10  10    % x
       2  1  0   2    1    0   2  1  0   2    1    0    2   1   0 ]; % y


A_e=inline('.5*((xe2*ye3-xe3*ye2)-(xe1*ye3-xe3*ye1)+(xe1*ye2-xe2*ye1))'...
    ,'xe1','xe2','xe3','ye1','ye2','ye3');

Be_11=inline('(ye2-ye3)','ye2','ye3');
Be_12=inline('(ye3-ye1)','ye1','ye3');
Be_13=inline('(ye1-ye2)','ye1','ye2');
Be_21=inline('(xe3-xe2)','xe2','xe3');
Be_22=inline('(xe1-xe3)','xe1','xe3');
Be_23=inline('(xe2-xe1)','xe1','xe2');

elx = zeros(n_el,4);                         % element x-coords var
ely = zeros(n_el,4);                         % element y-coords var

K = zeros(dof , dof);
f_Omega = zeros(dof ,1);
f_Gamma = zeros(dof ,1);

psi1GP = [.1666666666 .6666666666 .1666666666];
psi2GP = [.1666666666 .1666666666 .6666666666];
WGP = [.1666666666 .1666666666 .1666666666];

for e = 1:n_el
    Kmu_e = zeros(6);
    KL_e = zeros(6);

    for i=1:3
        elx(e,i) = GA(1,ICA(e,i));
        ely(e,i) = GA(2,ICA(e,i));
        
        if i==1
            elx(e,4) = GA(1,ICA(e,i));
            ely(e,4) = GA(2,ICA(e,i));
        end
    end
    xy_e = [elx(e,1:3)',ely(e,1:3)'];
    for j = 1:3
        psi1 = psi1GP(j);
        psi2 = psi2GP(j);
        psi3 = 1-psi1-psi2;
        
        N1 = psi1;
        N2 = psi2;
        N3 = psi3;  

        GN_e = [1 0 -1
                0 1 -1];

        Je = GN_e*xy_e;

        BB = Je\GN_e;

        Be = [BB(1,1),0,BB(1,2),0,BB(1,3),0;
               0,BB(2,1),0,BB(2,2),0,BB(2,3);
               BB(2,1),BB(1,1),BB(2,2),BB(1,2),BB(2,3),BB(1,3)];

        Le = [BB(1,1) BB(2,1), BB(1,2) BB(2,2), BB(1,3) BB(2,3)];           

        Kmu_e = Kmu_e + WGP(j)*(Be')*Dmu*Be*det(Je);
        KL_e = KL_e + WGP(j)*lamb*(Le)'*Le*det(Je);
        hold on
    end
        K(dofICA(e,1:6),dofICA(e,1:6)) = K(dofICA(e,1:6),dofICA(e,1:6)) + ...
          Kmu_e([1 3 5 2 4 6],[1 3 5 2 4 6]) + KL_e([1 3 5 2 4 6],[1 3 5 2 4 6]);
end

f_Gamma(25) = -f_*1/3;
f_Gamma(29) = f_*1/3;

f = f_Gamma;

K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;

K(3,:) = 0;
K(:,3) = 0;
K(3,3) = 1;

K(5,:) = 0;
K(:,5) = 0;
K(5,5) = 1;

K(6,:) = 0;
K(:,6) = 0;
K(6,6) = 1;
% Calculate unknown nodal tempretures using partitioning method
for el = 1:n_el

% patch(elx(e,1:3),ely(e,1:3),d(ICA(e,1:3)))
    for i=1:3
        elx(el,i) = elx(el,i);
        ely(el,i) = ely(el,i);

        if i==1
            elx(el,4) = elx(el,4);
            ely(el,4) = ely(el,4);
        end
    end
    plot (elx(el,1:4), ely(el,1:4),'Color','b')
    set(gcf, 'Visible', 'off')
    hold on
%     plot (elx(e,1:4), ely(e,1:4),'Color','w')
end
set(gcf, 'Visible', 'on')
d = K\f;

for el = 1:n_el

% patch(elx(e,1:3),ely(e,1:3),d(ICA(e,1:3)))
    for i=1:3
        elx(el,i) = elx(el,i)+ d(ICA(el,i)*2-1);
        ely(el,i) = ely(el,i) +d(ICA(el,i)*2);

        if i==1
            elx(el,4) = elx(el,4)+d(ICA(el,i)*2-1);
            ely(el,4) = ely(el,4)+d(ICA(el,i)*2);
        end
    end
    plot (elx(el,1:4), ely(el,1:4),'Color','r')
    set(gcf, 'Visible', 'off')
    hold on
%     plot (elx(e,1:4), ely(e,1:4),'Color','w')
end
set(gcf, 'Visible', 'on')
% title('Heat Conduction Triangular Element FEM'...
%             ,'FontWeight','normal','FontSize',14)
% % Create xlabel
% xlabel({'x (m)'},'FontWeight','normal','FontSize',14);
% % Create ylabel
% ylabel({'y (m)'},'FontWeight','normal'...
%     ,'FontSize',14);
% colorbar
% legend
x= 10;
y=2;
u_xy = 2*f_*(1-v^2)/(E*height)*x*(height/2-y);
disp('Analytical x displacement = ')
disp(u_xy)
disp('FE x displacement = ')
disp(d(13*2-1))
v_xy = f_*(1-v^2)/(E*height)*(x^2+v/(1-v)*y*(y-height));
disp('Analytical y displacement = ')
disp(v_xy)
disp('FE y displacement = ')
disp(d(13*2))