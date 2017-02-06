% Timo Comparison
%               e1      e2       e3
Exa  =      [0.02800 0.02658 0.02503];
Q1RI =      [0.02045 0.02018 0.01976];
Q1P0 =      [0.01798 0.01889 0.01974];
Q1   =      [0.01678 0.01373 7.699e-4];
P2   =      [0.02803 0.02661 0.02507];
P1P0 =      [3.160e-4 0.001990 0.004138];
P1   =      [.008793 0.007663 0.004317];
nu  =       [0.4        0.45    0.499];

% figure(1)
% plot(nu,Exa,nu,Q1RI,nu,Q1P0,nu,Q1)
% axis([0.4 0.5 0 0.03])
% xlabel({'\epsilon'},'FontWeight','demi','FontSize',14);
% % Create ylabel
% ylabel({'Displacement'},'FontWeight','demi'...
%     ,'FontSize',14);
% legend('Exact', 'Q1RI', 'Q1P0', 'Q1')
% print('Q', '-dpng', '-r600');
% 
% figure(2)
% plot(nu,Exa,nu,P2,nu,P1P0,nu,P1)
% axis([0.4 0.5 0 0.03])
% xlabel({'\gamma'},'FontWeight','demi','FontSize',14);
% % Create ylabel
% ylabel({'Displacement'},'FontWeight','demi'...
%     ,'FontSize',14);
% legend('Exact', 'P2', 'P1P0', 'P1')
% print('P', '-dpng', '-r600');

figure(3)
plot(nu,Exa,nu,Q1RI,nu,Q1P0,nu,Q1,nu,P2,nu,P1P0,nu,P1)
axis([0.4 0.5 0 0.03])
xlabel({'\nu'},'FontWeight','bold','FontSize',12);
% Create ylabel
ylabel({'y-displacement'},'FontWeight','demi'...
    ,'FontSize',12);
legend('Exact', 'Q1RI', 'Q1P0', 'Q1', 'P2', 'P1P0', 'P1')
print('All', '-dpng', '-r600');