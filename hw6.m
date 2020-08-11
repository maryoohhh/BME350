% BME 350 HW 6
% Mary Oh
% 1208315416

clc;
close all;
clear all;

%% Q6

% Transfer Function: H(s) = (s^2+10s+5)/(s^3+4s^2+10s+6)

% Numerator
num = [1 10 5];

% Denominator
den = [1 4 10 6];

% Tranfer function
H = tf(num, den)

%% Q6a

% Plot frequency response
figure(1)
bode(H), grid

%% Q6b

% Compute and plot poles and zero
[num, den] = eqtflength(num, den);
[z,p,k] = tf2zp(num, den)

figure(2)
pzplot(H), grid
text(real(z)+.1,imag(z),'Zero')
text(real(p)+.1,imag(p),'Pole')

%% Q6c

% Compute and plot step response
figure(3)
step(H)

%% Q6d

% Compute and plot impulse response
figure(4)
impulse(H)