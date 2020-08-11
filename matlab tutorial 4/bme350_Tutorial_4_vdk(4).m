
%%

%Exercise 1
clear all;
close all;
load kspace.mat; %Read in original kspace data
matx = 128;
maty = 128;

im = abs(fft2(kspace));
figure(1);imshow(abs(kspace), [0 2]);colormap jet;  %plot kspace

figure(2);imagesc(im);colormap gray; %plot image
%figure(20);surf(abs(kspace));%plot surface plot of kspace
%%

% Generate a 2D Fermi filter with radius "width" and transition region
% "edgewidth"
width = 30; % can be varied
edgewidth = 1; % can be varied
A = cumsum(ones(128,128),1);
B = cumsum(ones(128,128),2);
R = sqrt((matx/2.-A).*(matx/2.-A) + (maty/2.-B).*(maty/2.-B));
fermifilter =  1./(1+exp((R-width/2)/edgewidth));     


figure(3);surf(fermifilter);    %surface plot of fermifilter
figure(4);imshow(fermifilter, [0 1]);colormap jet; colorbar; %top view of fermifilter
%%

kspace_filtered = kspace.*fermifilter;% applying the filter
im2_filtered=abs(fft2(kspace_filtered));% filtered image

figure(5);imshow(abs(kspace_filtered), [0 2]);colormap jet;% display filtered kspace
figure(6);imshow(im2_filtered, [0 100]);colormap gray;% display filtered image


%%
%%Exercise 2

clear all
close all
clc

% Load image file
info = dicominfo('IM-0001-0001.dcm'); % Reads information of DICOM file. Change file name according to what you have.
imag1 = dicomread(info); % Get image information
imag1 = double(imag1); % Turn values from unsigned integers to double numbers 
% Display image
figure
imshow(abs(imag1));

 imcontrast % Allows you to dynamically change constrast of image
% imtool(imag1) % Allows you perform basic operations on the image.
%% FILTERING IMAGE

filter = 7; % Change variable to choice of 
hsize = 5; % Size of filter matrix

switch filter
    case 1
        h = fspecial('average',hsize); % Moving average filter
    case 2
        sigma = 5;
        h = fspecial('gaussian',hsize,sigma); % Low pass filter
    case 3
        sigma = 0.5;
        h = fspecial('log',hsize,sigma); % Edge detection
    case 4
        h = fspecial('prewitt'); % Vertical edge detection
    case 5
        h = fspecial('sobel'); % Horizontal edge detection.
    case 6
        radius = 3;
        h = fspecial('disk',radius); % Circular average filter
    case 7
        h = fspecial('unsharp',0.2); % Improves image resolution
end

imag2 = imfilter(imag1,h);
figure, imshow(imag2);imcontrast;
%% 
%Exercise 3
clc;
clear all;
close all;

% Load brain or knee data
load brain.mat
%load knee.mat


kspace = fftshift(fft2(im)); % Apply 2D fourier transform to obtain kspace data . 
figure(1);imshow(abs(im), [0 100]);colormap gray; % plot image. Change 100 to 5000 for knee data
colorbar;
title('Original Brain Image');

%%
% Exercise 3: Undersample kspace
kspace_us = kspace;   
kspace_us(2:2:end,:) = 0;  % Set every alternate row to zero

% Reconstruct image from undersampled kspace
im_us = ifft2(ifftshift(kspace_us));
figure(2);imshow(abs(im_us),[0 100]);colormap gray; % plot new image 
colorbar;
title('Brain Image from Undersampled kpace(every alternate line=0)');