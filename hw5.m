%% 4b

% C = 5uF
% R = 48 kOhm
% f = 6.63*10^-4 Hz

close all;
clear all;
clc;

% Transfer function: H(s)=0.240s/0.240s+1

% Numerator
num = [0.240];

% Denominator
den = [0.240 1];

% Transfer function
G = tf(num, den)

% Plot frequency response
figure (1)
bode(G), grid

%% 4c

% Transfer function: H(s)=0.240s/(0.53*10^-4s^2+0.242s+1)

% Numerator
num1 = [0.240];

% Denominator
den1 = [(0.53*10^-4) 0.242 1];

% Transfer function
G1 = tf(num1, den1)

% Plot frequency response
figure (2)
bode(G1), grid

%% 4d

% Transfer function: H(s)=(6.37*10^-4)s/(6.37*10^-4s+1)

% Numerator
num2 = [(6.37*10^-4)];

% Denominator
den2 = [(6.37*10^-4) 1];

% Transfer function
G2 = tf(num2, den2)

% Plot frequency response
figure (3)
bode(G2), grid
