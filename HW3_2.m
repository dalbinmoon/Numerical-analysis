%%
clear all ; close all ; clc;
%%
t = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
v  = [9.7 8.1 6.6 5.1 4.4 3.7 2.8 2.4 2.0 1.6 1.4 1.1 0.85 0.69 0.6];
figure()
plot(t, v); title('v verses t graph'); xlabel('t'); ylabel('v');
%%
R = 5000000;
H = [-2/R 1 ; -4/R 1 ; -6/R  1; -8/R 1 ; -10/R 1 ; -12/R 1 ; -14/R 1; -16/R 1 ; -18/R 1 ; -20/R 1 ; -22/R 1 ; -24/R 1 ; -26/R 1 ; -28/R 1 ; -30/R 1]
y  = [log(9.7) ; log(8.1) ; log(6.6) ; log(5.1) ; log(4.4) ; log(3.7) ; log(2.8) ; log(2.4) ; log(2.0) ; log(1.6) ; log(1.4) ; log(1.1) ; log(0.85) ; log(0.69) ; log(0.6)]

x = inv(transpose(H)*H)*(transpose(H)*y)
C = 1/x(1,1)
V = exp(x(2,1))
%%
x_line = linspace(2, 30);
y_line = V*exp(-x_line/(R*C));

figure()
line(x_line, y_line, 'Color', 'black'); hold on;
plot(t, v, 'o', 'Color', 'green');
title('LMS'); xlabel('t'); ylabel('v'); 
legend('origianl data', 'LMS line');