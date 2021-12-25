%% init
clear, clc
close all

I = imread('coins.png');

%% Edge Dection: Laplacian
% parameters setting
threshold = 80;
kernel = [-1, -1, -1
    -1, 8, -1;
    -1, -1, -1];

% perform filtering
E_La = conv2(I, kernel, 'same');
E_La = abs(E_La) > threshold;

% post-processing (optional)
se = strel('disk', 1);
E_La = imclose(E_La, se);
E_La = imopen(E_La, se);
E_La = imerode(E_La, se);

figure('Name', 'Laplacian'),
subplot(1,2,1),
imshow(I), title('Original')
subplot(1,2,2),
imshow(E_La), title('Edge')

%% Edge Detection: Sobel
% central difference = I(x+1, y) - I(x-1, y)
% this method is derived from central difference

threshold = 230;

% define kernels
% can also add diagonal kernals (optional)
% horizontal
k1 = [-1, 0, 1;
    -2, 0, 2;
    -1, 0, 1];
% vertical
k2 = rot90(k1);

% perform filtering
E1 = conv2(I, k1, 'same');
E2 = conv2(I, k2, 'same');
E_sobel = abs(E1) + abs(E2);
E_sobel = uint8(E_sobel);

figure,
subplot(1,2,1),
imshow(I), title('Original')
subplot(1,2,2),
imshow(E_sobel, [threshold 255]), title('Edge')

%% Edge Dection: Canny
% this method combine all steps
% post-processing: Linking
% strong(>threshold_high) + connected weak(>threshold_low) edge
E_canny = edge(I, 'canny');

figure,
subplot(1,2,1),
imshow(I), title('Original')
subplot(1,2,2),
imshow(E_canny), title('Edge')

%% Watershed (MATLAB func: watershed)
% concept: just like filling water into the container
% Pre-processing before Watershed:
% (1) distance transform: calculate distances from edge
% central point produce -> max
% points close to edge produce -> min
% (2) should perform "smoothing" before Watershed, reduce local min
% (3) marker-controlled segmentation
% step 1: marker image (self defined)
% define foreground and background
% step 2: minima imposition (MATLAB func: imimposemin)
% weaken local min (i.e. increase values), enhance global min (i.e. decrease values to 0)
% -----
% useful functions to design marker image
% imextendedmin
% imimcomplement

clear, clc
close all

rgb = imread('pears.png');
I = rgb2gray(rgb);
wat = watershed(I);

figure,
subplot(1,2,1),
imshow(I), title('Original')
subplot(1,2,2),
imshow(wat), title('Bad Watershed (w/o preprocessing)')

