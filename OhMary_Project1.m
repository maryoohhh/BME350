% Matlab Project 1
% Mary Oh
% 1208315416

%% 1.a Discrete time function
clear all;
clc;

asuid = 6;

% x[n]
n = -12:12;
x = zeros(1, length(n));

x(n >= -7 & n <= -1) = -2;
x(n >= 0 & n <= 6) =  0.5;

figure(1);
stem(n,x,'b'); %plot of discrete data
xlabel('n')
ylabel('x[n]')
xlim([-12 12])
ylim([-2.5 1])
title('Discrete time data - Mary Oh 1208315416')

%% 1.b Modified discrete time function

% x[n]
n = -20:20;
x = zeros(1, length(n));

x(n >= -7 & n <= -1) = -2;
x(n >= 0 & n <= 6) =  0.5;

figure(2);
stem(n,x,'b'); %plot of discrete data
xlabel('n')
ylabel('x[n]')
ylim([-3 2])
title('Modified discrete time data - Mary Oh 1208315416')

%% 1.c Reversed signal

% x[-n]
n = -20:20;
x = zeros(1, length(n));

x(n >= -7 & n <= -1) = -2;
x(n >= 0 & n <= 6) =  0.5;

figure(3);
stem(-n,x,'b'); %plot of discrete data
xlabel('-n')
ylabel('x[n]')
ylim([-3 2])
title('Reversed discrete time data - Mary Oh 1208315416')

%% 1.d Odd and even function

n = -20:20;
x = zeros(1, length(n));

x(n >= -7 & n <= -1) = -2;
x(n >= 0 & n <= 6) =  0.5;

xmt=[fliplr(x(n>=0)) fliplr(x(n<0))]

xe=0.5*(xmt+x)
xo=0.5*(x-xmt)

figure(4);
subplot(3,1,1);
stem(n,x);
xlabel('n')
ylabel('x[n]')
ylim([-3 2])
title('Signal - Mary Oh 1208315416')

subplot(3,1,2);
stem(n,xe);
xlabel('n')
ylabel('x[n]')
ylim([-3 2])
title('Even - Mary Oh 1208315416')

subplot(3,1,3);
stem(n,xo);
xlabel('n')
ylabel('x[n]')
ylim([-3 2])
title('Odd - Mary Oh 1208315416')

%% 2.a Continuous time signals

time = 0:0.01:10;

amp = 2;       % Amplitude [volts]
freq = 8.5;    % Frequency [radians/sec]
phase = -pi/3; % Phase     [radians]

% Define x
x = amp*cos(freq*time + phase); % Define x

% Plot
figure(5);
plot(time,x)
xlabel('Time (seconds)')
ylabel('Amplitude (volts)')
title('Continuous time signal, x(t) - Mary Oh 1208315416')
ylim([-2.5, 2.5]) % Adjusting the y-limits can make the plot easier to see

%% 2.b Period

% x(t) = 2*cos(8.5*t - pi/3). This signal is periodic.

freq = 8.5/(2*pi); % Hz
T = 1/freq

%% 2.c Sinusoidal signal y(t)

% Initialize

A = 6;
wo = ((asuid * 2) + 2)/(2*pi); % Hz
phi = pi/6;

T = 1/phi;
T1 = 18*T;

% Define y
time = 0:0.01:T1;
y = A*cos(wo*time+phi);

figure(6);
plot(time,y)

xlabel('Time (seconds)')
ylabel('Amplitude (volts)')
title('Sinusoidal time signal, x(t) - Mary Oh 1208315416')
ylim([-8, 8])

%% 3.a Discrete time signals

% Initialize
n = -20:20; % Time vector

x = zeros(1, length(n)); % Define x for every point in our time vector
h = zeros(1, length(n)); % Define h for every point in our time vector

% Assign values to x[n]
x(n >= -8 & n <= -6) = 1;
x(n >= -5 & n <= -3) = 0.5; 
x(n >= -2 & n <=  2) = 2;
x(n >=  3 & n <=  5) = 1;

figure(7);
stem(n,x,'k');
xlabel('n')
ylabel('x[n]')
ylim([-0.5 2.5])
title('Discrete time function - Mary Oh 1208315416')
hold on;

% create translated discrete step function
n0 = -(asuid+3);
n1 = 7;

u1 = n - n0;
u2 = n - n1;
%h = u1 - u2;

figure(8);
subplot(3,1,1);
stem(u1, x);
xlabel('n+asuid+3')
ylabel('u[n]')
ylim([-0.5 2.5])
title('Translated discrete step function - Mary Oh 1208315416')

subplot(3,1,2);
stem(u2, x, ':r');
xlabel('n-7')
ylabel('u[n]')
ylim([-0.5 2.5])
title('Translated discrete step function - Mary Oh 1208315416')

h1 = circshift(x, [n-n0]) - circshift(x, [n-n1]);
subplot(3,1,3);
stem(n, h1);
xlabel('n')
ylabel('h[n]')
ylim([-3 2.5])
title('Discrete impulse function - Mary Oh 1208315416')

%% Convolution function y[n]=x[n]*h[n]

y = conv(x,h1);

figure(9);
subplot(2,1,1);
stem(y);
xlabel('n')
ylabel('y[n]')
ylim([-25 25])
title('Convolution function x[n]*h[n] - Mary Oh 1208315416')

z = conv(h1,x)

figure(9);
subplot(2,1,2);
stem(z);
xlabel('n')
ylabel('y[n]')
ylim([-25 25])
title('Convolution function h[n]*x[n] - Mary Oh 1208315416')