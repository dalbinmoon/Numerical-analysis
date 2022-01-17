%%
close all; clear all; clc;
%% default value
P = 27.0;
D = 0.003;

zeta = 0.7;
%% acutator
omega_a = [2*pi*20.5 2*pi*16.4 2*pi*24.6];
num_a = [omega_a(2)^2];
den_a = [1 2*zeta*omega_a(2) omega_a(2)^2];
Ga = tf(num_a, den_a)
%% sensor
omega_s = 2*pi*100.0;
num_s = [omega_s^2];
den_s = [1 2*zeta*omega_s omega_s^2];
Gs = tf(num_s, den_s)
%% aerodynamics
Lp = 3.2;
Lpi = 1200;
Ldelta = 16000;
Gaero = tf(Ldelta, [1 Lp])
%% input output
theta_c = [pi/4 pi/8 0];
theta = theta_c()+1.0*(pi/180);
