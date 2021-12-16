%% init
clear, clc
close all

%%% parameters setting
% for image size (square)
N = 256;
% for unknown noise (Gaussian noise)
sigma = 0.05; % std for Gaussian
% for known noise (Motion noise)
len = 21;
theta = 13;
% for pseudo-inverse filtering
delta = 0.1;
% for limited inverse filtering
radius = N/2;

% read image
img = imresize(double(imread('cameraman.tif')), [N, N]);
img = (img-min(img(:)))/(max(img(:))-min(img(:)));
img = imnoise(img, 'gaussian', 0, sigma); % addicitve gaussian noise

% point spread function for motion noise
PSF = fspecial('motion', len, theta);

% generate motion-noised image
img_motion = imfilter(img, PSF, 'circular');

%% Section 1: Inverse filtering & Pseudo-inverse filtering & Radially-limited inverse filtering
% show images
figure,
subplot(3,2,1),
imshow(img), title('original image')
subplot(3,2,2),
imshow(img_motion), title('image with motion noise')

% frequency response of noise filter
PSF_fr = fft2(PSF, N, N);

% (a) inverse filtering
idx = abs(PSF_fr)>0; % prevent from 0 devided by 0 will produce NaN
img_invfilt = fft2(img_motion);
img_invfilt(idx) = img_invfilt(idx) ./ PSF_fr(idx);
img_invfilt = abs(ifft2(img_invfilt));
subplot(3,2,3),
imshow(img_invfilt), title('inverse filtered')

% (b) pseudo-inverse filtering
idx = abs(PSF_fr)>delta;
img_pseudo_invfilt = fft2(img_motion);
img_pseudo_invfilt(idx) = img_pseudo_invfilt(idx) ./ PSF_fr(idx);
img_pseudo_invfilt = abs(ifft2(img_pseudo_invfilt));
subplot(3,2,4),
imshow(img_pseudo_invfilt), title('pseudo-inverse filtered')

% (c) radially limited inverse filtering
idx = fftshift(circle_mask(img_motion, radius));
img_lim_invfilt = fft2(img_motion);
img_lim_invfilt(idx) = img_lim_invfilt(idx) ./ PSF_fr(idx);
img_lim_invfilt = abs(ifft2(img_lim_invfilt));
subplot(3,2,5), 
imshow(fftshift(idx)), title('radially limiting mask')
subplot(3,2,6)
imshow(img_lim_invfilt), title('radially limited inverse filtering')

%% Section 2: Deblur image using Wiener filter
% MATLAB function: "deconvwnr"
estimated_nsr = 0;
wnr2 = deconvwnr(img_motion, PSF, estimated_nsr);

figure,
subplot(1,2,1),
imshow(img_motion), title('image with motion noise')
subplot(1,2,2),
imshow(wnr2), title('restoration using deconvwnr')

%% Section 3: Deblur image using 
iter = 5;
luc1 = deconvlucy(img_motion, PSF, iter);

figure,
subplot(1,2,1),
imshow(img_motion), title('image with motion noise')
subplot(1,2,2),
imshow(wnr2), title('restoration using deconvlucy')

%% helper functions
function mask = circle_mask(img, radius)
% the processing region in pseudo-inverse filtering
mask = false(size(img));
c = round(size(mask)/2);
[x, y] = meshgrid(1:size(mask,2), 1:size(mask,1));
d = sqrt(((x-c(2)).^2+(y-c(1)).^2));
mask(d<=radius) = true;
end
