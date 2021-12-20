%% init
clear, clc
close all

%---parameters setting---
% for image size (256 by 256)
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
img = imnoise(img, 'gaussian', 0, sigma); % addicitve unknown gaussian noise

% point spread function for motion noise
PSF = fspecial('motion', len, theta);

% generate motion-noised image
img_motion = imfilter(img, PSF, 'circular');

%% Section 1: Inverse filtering & Pseudo-inverse filtering & Radially-limited inverse filtering
% show images
figure('Name', 'Section 1: Basic inverse filtering'),
subplot(3,2,1),
imshow(img), title('original image')
subplot(3,2,2),
imshow(img_motion), title('image with motion noise')

% frequency response of noise filter
PSF_fr = fft2(PSF, N, N);

% (a) inverse filtering
idx = abs(PSF_fr)>0; % prevent from 0 devided by 0 (this will produce NaN)
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

figure('Name', 'Section 2: Deblur by Wiener filter'),
subplot(1,2,1),
imshow(img_motion), title('image with motion noise')
subplot(1,2,2),
imshow(wnr2), title('deblured using deconvwnr')

%% Section 3: Deblur image using deconvlucy
% MATLAB function: "deconvlucy"
iter = 20;
lucy = deconvlucy(img_motion, PSF, iter);

figure('Name', 'Section 3: Deblur by deconvlucy'),
subplot(1,2,1),
imshow(img_motion), title('image with motion noise')
subplot(1,2,2),
imshow(lucy), title('deblured using deconvlucy')

%% Section 4: Deblur image using other methods
% MATLAB function: deconvblind
[blind, PSFr] = deconvblind(img_motion, PSF);

% MATLAB function: deconvreg
reg = deconvreg(img_motion, PSF);

figure('Name', 'Section 4: Deblur using other methods')
subplot(1,3,1),
imshow(img_motion), title('image with motion noise')
subplot(1,3,2),
imshow(blind), title('deblured using deconvblind')
subplot(1,3,3),
imshow(reg), title('deblured using deconvreg')

%% helper functions

function mask = circle_mask(img, radius)
% the processing region in pseudo-inverse filtering
mask = false(size(img));
c = round(size(mask)/2);
[x, y] = meshgrid(1:size(mask,2), 1:size(mask,1));
d = sqrt(((x-c(2)).^2+(y-c(1)).^2));
mask(d<=radius) = true;
end
