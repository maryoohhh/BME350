%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Project 4
%   Author: Mary Oh
%   Date: 2016/11/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Question 1a
clc;
clear all;
close all;
ASU = 6; % Input the last digit of your ASU ID number.

% Load brain data
load brain.mat
matx = 128;
maty = 128;

kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data

figure(1)
imshow(abs(im), [0 100]);
colormap jet; % plot image, adjust limits as needed
colorbar;
title('Original Brain Image')

figure(2)
imagesc(abs(kspace), [0 3.5e4]);
colormap jet;
colorbar;
title('Kspace Magnitude ASUID = 4');

figure(3) 
phasekspace = angle(kspace); 
imshow(phasekspace,[min(min(phasekspace)) max(max(phasekspace))]); 
title('Phase Brain Image ASUID = 6');

fprintf('The intensity range of the kspace is 3.5e4.\n')
fprintf('The coordinates of the point is about (64, 64).\n\n')

%% Question 1b

load brain.mat
matx = 128;
maty = 128;

kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data
kspace(64,76)=1e6;
brain=ifft2(kspace);

figure(4)
imshow(abs(brain), [0 100]);
colormap gray; % plot image
colorbar;
title('Original Brain Image 2');

figure(5)
imagesc(abs(kspace));
colormap jet;
title('Kspace Magnitude (64, 76) ASUID = 6');
colorbar;

fprintf('The image is lighter in color with streaks going downward from right to left.\n')
fprintf('The coordinates for this is about (64, 76).\n\n')

%%  Question 1c

load brain.mat
matx = 128;
maty = 128;

kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data
kspace(60,72)=1e6;
brain=ifft2(kspace);

figure(6)
imshow(abs(brain), [0 100]);
colormap gray; % plot image
colorbar;
title('Original Brain Image 3');

figure(7)
imagesc(abs(kspace));
colormap jet;
title('Kspace Magnitude (60,72) ASUID = 6');
colorbar;

fprintf('The image is lighter in color with streaks going downward from left to right.\n')
fprintf('The coordinates for this is about (60, 72).\n\n')

%%  Question 1d

load brain.mat
matx = 128;
maty = 128;

kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data
kspace(84,76)=1e6;
brain=ifft2(kspace);

figure(8)
imshow(abs(brain), [0 100]);
colormap gray; % plot image
colorbar;
title('Original Brain Image 4');

figure(9)
imagesc(abs(kspace));
colormap jet;
title('Kspace Magnitude (84,76) ASUID = 6');
colorbar;

fprintf('The image is lighter in color with streaks going downward from right to left.\n')
fprintf('The coordinates for this is about (84, 76).\n\n')

%% Question 2

clear all;
close all;
load kspace.mat; %Read in original kspace data

matx = 128;
maty = 128;
ASU = 6; % Input the last digit of your ASU ID number.

im = abs(fft2(kspace)); % Apply 2D fourier transform to obtain image data

figure(10)
imshow(abs(kspace), [0 2]);
colormap jet; % plot kspace
title('kspace');

figure(11)
imshow(im, [0 100]);
colormap gray; % plot image
title('Reconstructed Image');

%% Question 2a
% Generate a 2D Fermi filter with radius "width" and transition region
% "edgewidth"

width1 = (ASU+2)*2;
width2 = (ASU+2)*4;
width3 = (ASU+2)*8;
edgewidth = 4; % can be varied

A = cumsum(ones(128,128),1);
B = cumsum(ones(128,128),2);
R = sqrt((matx/2.-A).*(matx/2.-A) + (maty/2.-B).*(maty/2.-B));
fermifilter1 =  1./(1+exp((R-width1/2)/edgewidth));   
fermifilter2 =  1./(1+exp((R-width2/2)/edgewidth)); 
fermifilter3 =  1./(1+exp((R-width3/2)/edgewidth)); 

figure(12)
surf(fermifilter1); %surface plot of fermifilter

figure(13)
imshow(fermifilter1, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filtered1 = kspace.*fermifilter1;% applying the filter
im2_filtered1=abs(fft2(kspace_filtered1));% filtered image

figure(14)
imshow(abs(kspace_filtered1), [0 2]);
colormap jet;% display filtered kspace

figure(15)
imshow(im2_filtered1, [0 100]);
colormap gray;% display filtered image

figure(16)
surf(fermifilter2); %surface plot of fermifilter

figure(17)
imshow(fermifilter2, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filtered2 = kspace.*fermifilter2; % applying the filter
im2_filtered2=abs(fft2(kspace_filtered2));% filtered image

figure(18)
imshow(abs(kspace_filtered2), [0 2]);
colormap jet;% display filtered kspace

figure(19)
imshow(im2_filtered2, [0 100]);
colormap gray;% display filtered image

figure(20)
surf(fermifilter3); %surface plot of fermifilter

figure(21)
imshow(fermifilter3, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filtered3 = kspace.*fermifilter3;% applying the filter
im2_filtered3=abs(fft2(kspace_filtered3));% filtered image

figure(22)
imshow(abs(kspace_filtered3), [0 2]);
colormap jet;% display filtered kspace

figure(23)
imshow(im2_filtered3, [0 100]);
colormap gray;% display filtered image

%% Question 2b

kspace(84,76)=0.1;% Removing old spike. Change the coordinates according to location of spike
kspace(80+ASU,80+ASU)=50;% Adding new spike. Change the value according to intensity of spike

figure(24)
imshow(abs(kspace), [0 2]);
colormap jet; % plot kspace
title('kspace');

figure(25)
imshow(im, [0 100]);
colormap gray; % plot image
title('Reconstructed Image');

%% Question 2c

f1 = 2*(ASU+80-matx/2+16);
f2 = 2*(ASU+80-matx/2+2);
edgewidth = 1; % can be varied

A = cumsum(ones(128,128),1);
B = cumsum(ones(128,128),2);
R = sqrt((matx/2.-A).*(matx/2.-A) + (maty/2.-B).*(maty/2.-B));
fermifilterf1 =  1./(1+exp((R-f1/2)/edgewidth));   
fermifilterf2 =  1./(1+exp((R-f2/2)/edgewidth)); 

figure(26)
surf(fermifilterf1); %surface plot of fermifilter

figure(27)
imshow(fermifilterf1, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filteredf1 = kspace.*fermifilterf1;% applying the filter
im2_filteredf1=abs(fft2(kspace_filteredf1));% filtered image

figure(28)
imshow(abs(kspace_filteredf1), [0 2]);
colormap jet;% display filtered kspace

figure(29)
imshow(im2_filteredf1, [0 100]);
colormap gray;% display filtered image

figure(30)
surf(fermifilterf2); %surface plot of fermifilter

figure(31)
imshow(fermifilterf2, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filteredf2 = kspace.*fermifilterf2; % applying the filter
im2_filteredf2=abs(fft2(kspace_filteredf2));% filtered image

figure(32)
imshow(abs(kspace_filteredf2), [0 2]);
colormap jet;% display filtered kspace

figure(33)
imshow(im2_filteredf2, [0 100]);
colormap gray;% display filtered image

%% Question 2d

fermifilterf3 =  fermifilterf1 - fermifilterf2;   
fermifilterf4 =  fermifilterf1 - fermifilterf3;

figure(34)
surf(fermifilterf3); %surface plot of fermifilter

figure(35)
imshow(fermifilterf3, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filteredf3 = kspace.*fermifilterf3;% applying the filter
im2_filteredf3=abs(fft2(kspace_filteredf3));% filtered image

figure(36)
imshow(abs(kspace_filteredf3), [0 2]);
colormap jet;% display filtered kspace

figure(37)
imshow(im2_filteredf3, [0 100]);
colormap gray;% display filtered image

figure(38)
surf(fermifilterf4); %surface plot of fermifilter

figure(39)
imshow(fermifilterf4, [0 1]);
colormap jet;
colorbar; %top view of fermifilter

kspace_filteredf4 = kspace.*fermifilterf4; % applying the filter
im2_filteredf4=abs(fft2(kspace_filteredf4));% filtered image

figure(40)
imshow(abs(kspace_filteredf4), [0 2]);
colormap jet;% display filtered kspace

figure(41)
imshow(im2_filteredf4, [0 100]);
colormap gray;% display filtered image

%% Question 3
clear all;
close all;
clc;
ASU = 6; % Input the last digit of your ASU ID number.
% Load image file
info = dicominfo('IM-0001-0001.dcm'); % Reads information of DICOM file. Change file name according to what you have.
imag1 = dicomread(info); % Get image information
imag1 = double(imag1); % Turn values from unsigned integers to double numbers
imag1 = imag1/max(max(imag1)); 

% Display image
figure(42)
imshow(imag1);
imcontrast % Allows you to dynamically change constrast of image

%% Question 3a
filter = 4; % Change variable as per question to choose filter.
hsize = ASU+1; % Size of filter matrix. Change as per question.

switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = 8;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 0.4;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Horizontal edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end

imag2 = imfilter(imag1,h);
figure(43)
imshow(imag2);
title('1208315416, hsize 7');

filter = 4; % Change variable as per question to choose filter.
hsize = (ASU+1)*2; % Size of filter matrix. Change as per question.

switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = 8;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 0.4;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Horizontal edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end

imag2 = imfilter(imag1,h);
figure(44)
imshow(imag2);
title('1208315416, hsize 14');

filter = 4; % Change variable as per question to choose filter.
hsize = (ASU+1)*8; % Size of filter matrix. Change as per question.

switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = 8;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 0.4;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Horizontal edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end

imag2 = imfilter(imag1,h);
figure(45)
imshow(imag2);
title('1208315416, hsize 56');

%% Question 3b

info = dicominfo('IM-0001-0001.dcm'); % Reads information of DICOM file. Change file name according to what you have.
imag1 = dicomread(info); % Get image information
imag1 = double(imag1); % Turn values from unsigned integers to double numbers
imag1 = imag1/max(max(imag1)); 

filter = 4;
hsize = 8;

switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = (ASU+10)/2;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 8.5;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Horizontal edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end
imaga = imfilter(imag1,h);
figure(46), imshow(imaga);
title ('Gaussian Filtered Image');

filter = 3;
hsize = 8;
switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = 8.5;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 0.2;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Horizontal edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end
imagb = imfilter(imaga,h);
figure(47), imshow(imagb);
title ('Gaussian Filtered & Edge Detection Image');

%% Question 4
% Load and display data
clc;
clear all;
close all;
ASU = 6; % Input the last digit of your ASU ID number.

load brain.mat % Load brain  data
kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data . 

figure(48)
imshow(abs(im), [0 100]);
colormap gray; % plot image. 
colorbar;
title('Original Brain Image');

%% Question 4a

% Undersample kspace
kspace_us = kspace;   
kspace_us(2:2:end,:) = 0;  % Setting every alternate row to zero. 

figure(49)
subplot(3,1,1)
imshow(abs(kspace), [0 100]);
colormap gray; % plot image. 
colorbar;
title('Original kspace');

subplot(3,1,2)
imshow(abs(im), [0 100]);
colormap gray; % plot image. 
colorbar;
title('Original Brain Image');

% Reconstruct image from undersampled kspace
im_us = ifft2(ifftshift(kspace_us));
subplot(3,1,3)
imshow(abs(im_us),[0 100]);
colormap gray; % plot new image. You may have to decrease the upper limit, i.e. 100 , to appreciate what is going on. 
colorbar;
title('Reconstructed Brain Image from Undersampled kpace(every alternate line=0)');

%% Question 4b
