%% discrete sin functions

n = 10; % changes fraction pi/n 
t = 0:1/n:10; %counting by pi/n

figure;

ax1 = subplot(4,1,1);
w0 = 1; %w0 = angular frequency terms of 2pi;
y = cos(t*pi/w0);
stem(t,y,'k');

ax2 = subplot(4,1,2);
w0 = 2;
y = cos(t*pi/w0);
stem(t,y,'b');

ax3 = subplot(4,1,3);
w0 = 4;
y = cos(t*pi/w0);
stem(t,y,'r');

ax4 = subplot(4,1,4);
w0 = 8;
y = cos(t*pi/w0);
stem(t,y,'g');

linkaxes([ax1,ax2,ax3,ax4],'x');


%% Odd and Even functions
% x[-n] = x[n] is Even
% x[-n] - x[n] = 0, test for "evenness"

% x[-n] = -x[n] is Odd
% x[-n] - -x[n] = 0, test for "oddness"

% creat triangle function
n1 = 0:pi/10:2*pi;
x1 = sawtooth(n1,0.5)+1;
n2 = -2*pi:pi/10:0;
x2 = sawtooth(n2,0.5)+1;
n = [n2 n1]; % concatenate n1 and n1 into n
x = [-x2 x1]; % concatenate x1 and x2 into x
plot(n,x,'k','LineWidth',2)
grid on;

% is x Even?
figure;
plot(n,x,'k','LineWidth',2)
hold on;
plot(-n,x,'.r','LineWidth',2,'MarkerSize',20)

% is x Odd?
figure;
plot(n,-x,'k','LineWidth',2)
hold on;
plot(-n,x,'.r','LineWidth',2,'MarkerSize',20)


%% Simple Decomposition

% x[n]
n = -5:5;
x = [0 0 0 0 0 1 1 1 1 1 1];
figure;
stem(n,x,'b'); %lollipop plot of discrete data

% x_even = (x[n] + x[-n])/2
%fliplf just flips the vector around the x-axis, same as x[-n]
figure;
x_even = (x + fliplr(x))/2;
stem(n,x_even,'k');

% x_odd = (x[n] - x[-n])/2
%fliplf just flips the vector around the x-axis, same as x[-n]
figure;
x_odd = (x - fliplr(x))/2;
stem(n,x_odd,'k');

% recompose (sanity check)
figure;
subplot(3,1,1);
stem(n,x_even,'b');

subplot(3,1,2);
stem(n,x_odd,'r');

% recreate original signal by adding x Even and Odd
x_recompose = x_even + x_odd; 
subplot(3,1,3);
stem(n,x_recompose,'k');


%% Energey and Power
clear;
clc;

% Energy
% function x[n]
amp = 10;
time = 100;
n = -time*pi:pi/100:time*pi;

x = amp * sin(n);

% ploting x[n]
figure;
plot(n,x,'k','LineWidth',2);
hold on;

% retify signal
plot(n,abs(x),':r','LineWidth',2);
% square rectified signal
plot(n,abs(x).^2,'b','LineWidth',2);

% calculate Energy
energy = sum(abs(x).^2) % tends to infinity for increasing n


% Power
% function x[n]
amp = 1;
time = 100;
n = -time*pi:pi/100:time*pi;

x = amp * sin(n);

% ploting x[n]
figure;
plot(n,x,'k','LineWidth',2);
hold on;

% retify signal
plot(n,abs(x),':r','LineWidth',2);
% square rectified signal
plot(n,abs(x).^2,':b','LineWidth',2);

% calculate Power
power = sum((abs(x).^2)/(length(n)+1)) % proportioanl to amp


%% Discrete Impulse function
clear;
clc;

% create discrete impulse function
n = -10:10; % generate time base
x = zeros(1,length(n)); % generate vector of zeros, same size at n
x(n==0) = 1; % set 0 time point = 1, to create impulse function

figure;
stem(n,x,'k','LineWidth',2);

sum(x(find(n==-10):find(n==-1))) % should = 0
sum(x(find(n==1):find(n==10))) % should = 0
sum(x) % should = 1

% moving impulse in time
% just like translation of any function
n0 = 4;
n_trans = n - n0;

figure;
stem(n,x,':b','LineWidth',2);
hold on;
stem(n_trans,x,':r','LineWidth',2);

% using impulse function to sample another function
x_sin = sin(n); % create sin function base on n
figure;
stem(n,x_sin,'k','LineWidth',2);
hold on;
stem(n-4,x,':b','LineWidth',2);
stem(n-3,x,':r','LineWidth',2);

% multiply sin fucntion by impulse function to sample
figure;
x_sample = x_sin .* circshift(x,[0 -4]); % shift impulse function by n0
stem(n,x_sample,'b','LineWidth',2);
hold on;
x_sample = x_sin .* circshift(x,[0 -3]); % shift impulse function by n0
stem(n,x_sample,'r','LineWidth',2);

%% Discrete Step function
clear;
clc;

% create discrete step function
n = -10:10; % generate time base
u = zeros(1,length(n)); % generate vector of zeros, same size at n
u(find(n==0):end) = 1; % set time point >= 0 to 1

figure;
stem(n,u,'k','LineWidth',2);
hold on;
% create translated discrete step function
n0 = -5;
stem(n - n0,u,':r','LineWidth',2);

% substract step functions to create impulse function
% d[n] = u[n] - u[n-n0] where n0 = 1
n0 = 1;
d = u - circshift(u,[0 n0]); % shift impulse function by n0
figure;
stem(n,d,'k','LineWidth',2);
% why is there a -1 impulse at -10???
% difference between math & code.

