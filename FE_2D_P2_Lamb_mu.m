clear all ; close all ; clc ;
npoints = 100;
% Set the number of elements
n_el = 16; n_np = 45;
dof = 2 * n_np;
f_ = 10;
height = 2;
Length = 10;              % length of bar
% Set the properties of the material and problem
E = 3000;          % Youngs Mod
v = .5;
% 
f = (1-v);    g = v;   a=E/((1+v)*(1-2*v)); q=(1-2*v)/2;
D = a*[f g 0;g f 0;0 0 q];
w = 6;          % weight of element

ICA = [1 3 13 2 8 7 % Interconnectivity array
       1 13 11 7 12 6
       3 5 15 4 10 9
       3 15 13 9 14 8
       11 13 23 12 18 17
       11 23 21 17 22 16
       13 15 25 14 20 19
       13 25 23 19 24 18
       21 23 33 22 28 27
       21 33 31 27 32 26
       23 25 35 24 30 29
       23 35 33 29 34 28
       31 33 43 32 38 37
       31 43 41 37 42 36
       33 35 45 34 40 39
       33 45 43 39 44 38];

dofICA = [ICA*2-1,ICA*2];
   

GA1=[0  0   0 0  0 1.25 1.25 1.25 1.25 1.25 2.5 2.5 2.5 2.5 2.5 3.75 3.75 3.75;
     2  1.5 1 0.5 0 2    1.5  1.   0.5  0    2   1.5 1   0.5 0    2   1.5 1.];
 
GA2 =[3.7500  3.7500  5.0000  5.0000  5.0000  5.0000  5.0000  6.2500  6.2500;
      0.5000    0     2.0000  1.5000  1.0000  0.5000     0    2.0000  1.5000];

GA3 =[6.2500  6.2500  6.2500  7.5000  7.5000  7.5000  7.5000  7.5000  8.7500;
      1.0000  0.5000   0      2.0000  1.5000  1.0000  0.5000     0    2.0000];

GA4 =[8.7500  8.7500  8.7500  8.7500  10.0000  10.0000  10.0000  10.0000  10;
      1.5000  1.0000  0.5000    0     2.0000   1.5000   1.0000   0.5000   0];

GA = [GA1, GA2, GA3, GA4];

A_e=inline('.5*((xe2*ye3-xe3*ye2)-(xe1*ye3-xe3*ye1)+(xe1*ye2-xe2*ye1))'...
    ,'xe1','xe2','xe3','ye1','ye2','ye3');

psi1GP = [.1666666666 .6666666666 .1666666666];
psi2GP = [.1666666666 .1666666666 .6666666666];
WGP = [.1666666666 .1666666666 .1666666666];

elx = zeros(n_el,7);                         % element x-coords var
ely = zeros(n_el,7);                         % element y-coords var

draw = [1 4 2 5 3 6 7 1];

Ke = zeros(12 , 12);
K = zeros(dof , dof);
f_Omega = zeros(dof ,1);
f_Gamma = zeros(dof ,1);

for e = 1:n_el
    Ke = zeros(12 , 12);
    for gpn = 1:3
        psi1 = psi1GP(gpn);
        psi2 = psi2GP(gpn);
        psi3 = 1-psi1-psi2;

        for i=1:6
            elx(e,i) = GA(1,ICA(e,i));
            ely(e,i) = GA(2,ICA(e,i));

            if i==1
                elx(e,7) = GA(1,ICA(e,i));
                ely(e,7) = GA(2,ICA(e,i));
            end
        end
        xy_e = [elx(e,1:6)',ely(e,1:6)'];

        N1 = psi1*(2*psi1-1);
        N2 = psi2*(2*psi2-1);
        N3 = psi3*(2*psi3-1);
        N4 = 4*psi1*psi2;
        N5 = 4*psi2*psi3;
        N6 = 4*psi1*psi3;               % Table 7.5 Pg 175 Belytschko
%%%%%%%%%%
        GN_e = [4*psi1-1, 0, -3+4*psi1+4*psi2, 4*psi2, -4*psi2, 4-4*psi2-8*psi1;
                0, 4*psi2-1, -3+4*psi2+4*psi1, 4*psi1, 4-4*psi1-8*psi2, -4*psi1];
%%%%%%%%%%%check
        Je = GN_e*xy_e;

        BB = Je\GN_e;
        Be = [BB(1,1),0,BB(1,2),0,BB(1,3),0,BB(1,4),0,BB(1,5),0,BB(1,6),0;
              0,BB(2,1),0,BB(2,2),0,BB(2,3),0,BB(2,4),0,BB(2,5),0,BB(2,6);
              BB(2,1),BB(1,1),BB(2,2),BB(1,2),BB(2,3),BB(1,3),BB(2,4),BB(1,4),...
              BB(2,5),BB(1,5),BB(2,6),BB(1,6)];

        Ke = Ke + WGP(gpn)*(Be')*D*Be*det(Je);
    end
        
    K(dofICA(e,1:12),dofICA(e,1:12)) = K(dofICA(e,1:12),dofICA(e,1:12)) + ...
        Ke([1 3 5 7 9 11 2 4 6 8 10 12],[1 3 5 7 9 11 2 4 6 8 10 12]);

    hold on

end

f_Gamma(41*2-1) = -f_*.2083;
f_Gamma(42*2-1) = -f_*.25;
f_Gamma(44*2-1) = f_*.25;
f_Gamma(45*2-1) = f_*.2083;

F = f_Gamma;

K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;

K(3,:) = 0;
K(:,3) = 0;
K(3,3) = 1;

K(5,:) = 0;
K(:,5) = 0;
K(5,5) = 1;

K(7,:) = 0;
K(:,7) = 0;
K(7,7) = 1;

K(9,:) = 0;
K(:,9) = 0;
K(9,9) = 1;

K(10,:) = 0;
K(:,10) = 0;
K(10,10) = 1;

% Calculate unknown nodal tempretures using partitioning method
for el = 1:n_el
    plot(elx(el,draw), ely(el,draw),'Color','b')
    set(gcf, 'Visible', 'off')
    hold on
end
set(gcf, 'Visible', 'on')
d = K\F;

for el = 1:n_el
    for i=1:6
        elx(el,i) = elx(el,i)+d(ICA(el,i)*2-1);
        ely(el,i) = ely(el,i)+d(ICA(el,i)*2);

        if i==1
            elx(el,7) = elx(el,7)+d(ICA(el,i)*2-1);
            ely(el,7) = ely(el,7)+d(ICA(el,i)*2);
        end
    end
    plot(elx(el,draw), ely(el,draw),'Color','r')
    set(gcf, 'Visible', 'off')
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
format long
u_xy = 2*f_*(1-v^2)/(E*height)*x*(height/2-y);
disp('Analytical x displacement = ')
disp(u_xy)
disp('FE x displacement = ')
disp(d(41*2-1))
v_xy = f_*(1-v^2)/(E*height)*(x^2+v/(1-v)*y*(y-height));
disp('Analytical y displacement = ')
disp(v_xy)
disp('FE y displacement = ')
disp(d(41*2))