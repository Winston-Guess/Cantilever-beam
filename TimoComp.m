% Timo Comparison
%               e1      e2       e3
Quad10 =    [0.01510 0.003854 0.002729];
Lin10 =     [0.01498 0.003558 0.001489];
Lin30 =     [0.01509 0.003819 0.002498];
Lin100 =    [0.01510 0.003851 0.002707];
LinUI10 =   [0.01500 0.003750 0.002625];
LinUI30 =   [0.01509 0.003843 0.002718];
LinUI100 =  [0.01510 0.003853 0.002728];
eps = [0.1 0.01 0.001];

figure(1)
loglog(eps,Quad10,eps,Lin10,eps,Lin30,eps,Lin100)
axis([0 0.1 0 0.02])
xlabel({'\epsilon'},'FontWeight','demi','FontSize',14);
% Create ylabel
ylabel({'Displacement'},'FontWeight','demi'...
    ,'FontSize',14);
legend('10 Quadratic Elements', '10 Linear Elements', '30 Linear Elements', '100 Linear Elements')
print('QL', '-dpng', '-r600');

figure(2)
loglog(eps,Quad10,eps,LinUI10,eps,LinUI30,eps,LinUI100)
axis([0 0.1 0 0.02])
xlabel({'\epsilon'},'FontWeight','demi','FontSize',14);
% Create ylabel
ylabel({'Displacement'},'FontWeight','demi'...
    ,'FontSize',14);
legend('10 Quadratic Elements', '10 Linear Elements', '30 Linear Elements', '100 Linear Elements')
print('QLUI', '-dpng', '-r600');