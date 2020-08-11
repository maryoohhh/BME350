%% simple function

n = -5:5;
x = [0 0 0 0 -1 1 1 2 1 0 0]; % think of these as the weights on detla
figure;
stem(n,x,'r','LineWidth',3); 

%% decomposing simple funtion

delta = [1]; %extremely simple delta function

%this is just vector of color options for ploting differnt colors
c = ['r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c'];

figure;
% decompose the funtion into separate translated and 
% weighted delta functions
for j = 1:length(x)
    y(j) = x(j) * delta
    
    subplot(3,4,j);
    stem(n(j),y(j),'LineWidth',2,'Color',c(j));
    ylim([-1 4]);
end;

% recompose the funtion from separate translated and 
% weighted delta functions
for j = 1:length(x)
    y(j) = x(j) * delta
    
    subplot(3,4,12);
    hold on
    stem(n(j),y(j),'LineWidth',2,'Color',c(j));
    ylim([-1 4]);
end;

%% simple convolution sum

% system
x = [0 0 0 0 0 1 1 1 1 1 3 3 3 3 3 -1 -1 -1 -1 -1 0 0 0 0 0]; 
% x = [0 0 0 0 0 1 1 1 1 1 0 0 0 0 0];

% impulse
% h = [1]; % simple edge detection function/matrix, try different values
h = [-1 1]; % simple edge detection function/matrix, try different values
% h = [1 0 -1]; % simple edge detection function/matrix, try different values
% h = [1 0.5 -0.5 -1]; % simple edge detection function/matrix, try different values
% h = [2 1 0 -1 -2]; 

h_flip = fliplr(h); % h = h(-t);

figure;
for n = 1:length(x)-(length(h)-1)
     
    subplot(3,1,1)
    cla; % clear current axes
    stem(x,'b','LineWidth',2); % plot x
    hold on;
    stem(n:n+length(h_flip)-1,h_flip,'r','LineWidth',2); % plot h
    xlim([0 20]);
    
    % plot y, "pseudo" convolution sum (doesn't capture leading and
    % trailing zeros)
    subplot(3,1,2)
    stem(n+length(h_flip)-1,sum(x(n:n+length(h_flip)-1) .* h_flip),'k','LineWidth',3);
    hold on;
    xlim([0 20]);
    
    pause;
end

% check with Matlab convoltion function
subplot(3,1,3)
y = conv(h,x); 
stem(y,'k','LineWidth',3);
xlim([0 20]);

%% 2-D convolution

A = zeros(10);
A(3:7,3:7) = ones(5);
figure;
mesh(A);

s = [1 2 1; 0 0 0; -1 -2 -1]; % edge detection function/matrix
% why does this fucntion find edges when convolved with another function?
figure;
mesh(s);

figure;
H = conv2(A,s); % convolve to detect horizontal edges
subplot(1,3,1);
mesh(H);
V = conv2(A,s'); % transpose and convolve to detect vertical edges
subplot(1,3,2);
mesh(V);
subplot(1,3,3);
mesh(sqrt(H.^2 + V.^2)); % combine and normalize results

%% convolution for image processing
% convolution forms basis for more advance methods
I = imread('circuit.tif');
figure;
imshow(I);
BW1 = edge(I,'canny');
figure;
imshow(BW1);


