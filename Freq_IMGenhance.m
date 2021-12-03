% ideal filter
% filter design, why we design the non-ideal filters?
% because: ideal filter in Freq Domain is NOT ideal in Spatial Domain, so
% we are doing the convolution of rippling-filter to my images. 

%% Ideal Filter
clear, clc
close all

sz_r = 257; % odd number, to produce the symmerticity
sz_c = 257; % odd number, to produce the symmerticity
radius = 10;
img = imresize(double((imread('moon.tif'))),[sz_r, sz_c]);

% design ideal low-pass filter
filter = zeros(sz_r, sz_c);
c = round(size(filter)/2);% center of the circle
[X, Y] = meshgrid(1:sz_c, 1:sz_r);
d = sqrt((X-c(2)).^2+(Y-c(1)).^2);
filter(d<radius') = 1;

% perfrom filtering
fimg = fftshift(fft2(img));
img_filt = ifft2(ifftshift(fimg.*filter));

% display
figure,
subplot(1,3,1), imagesc(img), title('original image')
subplot(1,3,2), imagesc(filter), title('filter, in frequency domain')
subplot(1,3,3), imagesc(img_filt), title('filtered image (ideal filter)')
colormap gray

%% Unsharp filtering
% unsharp filtering: f_sharp = f + w*(f - f_low-pass)
% functions:
% (1) fspecial: to create the filter
% (2) imfilter: to filter the image
clear, clc
close all

w = 1; % weight, for unsharp filtering
img = imread('tire.tif');

% low-pass filtering (using Gaussian filter)
filter = fspecial('average', 9); % a low-pass filter
img_lp = imfilter(img, filter);

% display
img_sharp = img + w*(img - img_lp);
img_sharp = imhistmatch(img_sharp, img);

figure,
subplot(1,3,1), imshow(img), title('original image')
subplot(1,3,2), imshow(img_lp), title('low-pass filtered image')
subplot(1,3,3), imshow(img_sharp), title('sharped image')

%% High boost filtering
% high boost filtering: f_sharp = (A-1)*f + f_high-pass
% functions:
% (1) fspecial: to create the filter
% (2) imfilter: to filter the image
clear, clc
close all

alpha = 0; % alpha for laplacian filter
A = 1.5; % weight, for high boost filtering
img = imread('tire.tif');

% low-pass filtering (using Gaussian filter)
filter = fspecial('laplacian', alpha); % a high-pass filter
img_hp = imfilter(img, filter);

% display
img_sharp = (A-1)*img + img_hp;
img_sharp = imhistmatch(img_sharp, img);

figure,
subplot(1,3,1), imshow(img), title('original image')
subplot(1,3,2), imshow(img_hp), title('high-pass filtered image')
subplot(1,3,3), imshow(img_sharp), title('sharped image')
