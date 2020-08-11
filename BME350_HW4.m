% BME 350 HW 4
% Mary Oh
% 1208315416

% x(t) = cos(30*pi*t)
% y(t) = sin(56*pi*t)
% z(t) = cos(72*pi*t)

Fs = 100;
Ss = 1/Fs;
t = -0.5:Ss:0.5

% original function
x = cos(30*pi*t);
y = sin(56*pi*t);
z = cos(72*pi*t);

x_NFFT = length(x);
y_NFFT = length(y);
z_NFFT = length(z);

x_freqbase = Fs/2*linspace(0,1,x_NFFT/2+1);
y_freqbase = Fs/2*linspace(0,1,y_NFFT/2+1);
z_freqbase = Fs/2*linspace(0,1,z_NFFT/2+1);

X = fft(x, x_NFFT);
Y = fft(y, y_NFFT);
Z = fft(z, z_NFFT);

figure(1)
subplot(3,1,1)
plot(x_freqbase, abs(X(1:x_NFFT/2+1)),'k')
title('x(t) = cos(30*pi*t)')
xlabel('Frequency (Hz)')
ylabel('Power')

subplot(3,1,2)
plot(y_freqbase, abs(Y(1:y_NFFT/2+1)),'k')
title('y(t) = sin(56*pi*t)')
xlabel('Frequency (Hz)')
ylabel('Power')

subplot(3,1,3)
plot(z_freqbase, abs(Z(1:z_NFFT/2+1)),'k')
title('z(t) = cos(72*pi*t)')
xlabel('Frequency (Hz)')
ylabel('Power')

%% part b

v = 2*x + 4*y - 3*z;

v_NFFT = length(v);
V = fft(v, v_NFFT);

v_freqbase = Fs/2*linspace(0,1,v_NFFT/2+1);

figure(2)
subplot(2,1,1)
plot(t,v)
title('v(t) = 2x(t) + 4y(t) - 3z(t)')
xlabel('Frequency (Hz)')
ylabel('Power')

subplot(2,1,2)
plot(v_freqbase, abs(V(1:v_NFFT/2+1)),'k')
title('Fourier transform')
xlabel('Frequency (Hz)')
ylabel('Power')