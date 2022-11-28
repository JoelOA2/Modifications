clear;
% Let:
% x(1) = Target cell (CD4+ T-cell population), T
% x(2) = Number of Infected Cells, I
% X(3) = Virus Concentration, V

% Independent variables:
b = 0.17;
dt = 0.01;
B = 6.5*10^-4;
di = 0.39;
c = 3;

% The differential equations:
f = @(t, x) [b - dt*x(1) - B*x(3)*x(1);
    B*x(3)*x(1) - di*x(2);
    x(4)*x(2) - c*x(3);
    -0.0077*x(4)];

T0 = 10; % Initial value of x(1) / T
I0 = 0;  % Initial value of x(2) / I
V0 = 10^-6;  % Initial value of x(3) / V
po = 850;

[t, xsol] = ode45(f, [0 500], [T0 I0 V0 po]);

figure();
yyaxis left;
semilogy(t, xsol(:, 3)*10^3 , 'r');
axis([0 500 1 10^8]);
ylabel('Virus Load');
title('Model Dynamics of CD4+ T-cell Count Against Virus Load');
xlabel('Time / (days)');
hold on

yyaxis right;
plot(t ,xsol(:, 1), 'b')
ylabel('T-Lymphocytes CD4+ Count');


ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';