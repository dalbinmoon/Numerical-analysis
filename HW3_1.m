%%
clear all ; close all ; clc;
%%
'1_(a)'
T = [0 10 20 30 40 50 60 70 80 90 100];
P = [0.94 0.96 1.00 1.05 1.07 1.09 1.14 1.17 1.21 1.24 1.28];
figure()
plot(T, P); title('P verses T graph'); xlabel('T'); ylabel('P');
%%
'1_(b) : Calculating LSM by hand'
H = [0 1 ; 30 1 ; 70 1; 100 1];
Y = [0.94; 1.05; 1.17; 1.28];

x = inv(transpose(H)*H)*(transpose(H)*Y)
%%
'1_(c)'
%LSM으로 구한 값으로 직선의 방정식 표현하기
x_line = linspace(0, 100);
y_line = x(1,1)*x_line + x(2,1);

%extrapolate 구하기
x_over = 100 : 5 : 150;
y_over = x(1,1)*x_over + x(2,1)

figure()
line(x_line, y_line, 'Color', 'black'); hold on;
plot(T, P,'o', 'Color', 'green'); hold on ; 
plot(x_over, y_over, '*', 'Color', 'red'); xlim([0 150]);
title('extrapolate'); xlabel('T'); ylabel('P'); 
legend('original line', 'origianl data', 'extrapolate data');
%%
'1_(c) : Finding To'
x_To = linspace(-500, 500);
y_To = x(1,1)*x_To + x(2,1);

To = -x(2,1)/x(1,1);

figure()
line(x_To, y_To, 'Color', 'black'); hold on;
plot(To, 0, '+', 'Color', 'red', 'LineWidth', 3);
title('Finding To'); xlabel('T'); ylabel('P'); 