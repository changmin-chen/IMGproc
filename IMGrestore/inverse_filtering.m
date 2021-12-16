%% init
clear, clc
close all

% parameters
delta = 0.1;
sigma = 0.01;
radius = 30;

% load image
img = double(imread('cameraman.tif'));
img = (img-min(img(:)))/(max(img(:))-min(img(:)));
img = imnoise(img, 'gaussian', 0, sigma); % unknown noise

% circle mask
mask = false(size(img));
c = round(size(mask)/2);
[x, y] = meshgrid(1:size(mask,2), 1:size(mask,1));
d = sqrt(((x-c(2)).^2+(y-c(1)).^2));
mask(d<=radius) = true;
imshow(mask)

%% Section 1: Inverse filtering & Pseudo-inverse filtering
% define kernel which represents noise in frequency domain
kernel = fspecial('gaussian', 256, 25);
kernel = mat2gray(kernel);

% apply kernel, generate noisy image
img_freq = fftshift(fft2(img));
img_freq_noise = img_freq.*kernel;
img_noise = abs(ifft2(fftshift(img_freq_noise)));

% try inverse filtering
img_invfilt = img_freq_noise ./ kernel;
img_invfilt(isnan(img_invfilt)) = 0;
img_invfilt = abs(ifft2(fftshift(img_invfilt)));

% try pseudo-inverse filtering
H = 1 ./ kernel;
idx = kernel > delta;
img_pinvfilt = zeros(size(img));
img_pinvfilt(idx) = img_freq_noise(idx).*H(idx);
img_pinvfilt = abs(ifft2(fftshift(img_pinvfilt)));

% try limited inverse filtering
mask_pseudo = mask & idx;
img_liminvfilt = zeros(size(img));
img_liminvfilt(mask_pseudo) = img_freq_noise(mask_pseudo).*H(mask_pseudo);
img_liminvfilt(isnan(img_liminvfilt)) = 0;
img_liminvfilt = abs(ifft2(fftshift(img_liminvfilt)));

% try "deconvwnr"

% display
figure,
subplot(2,3,1), imshow(img_noise), title('noisy image')
subplot(2,3,2), imshow(img_invfilt), title('inv filt')
subplot(2,3,3), imshow(img_pinvfilt), title('pseudo inv filt')
subplot(2,3,4), imshow(mask_pseudo), title('mask & pseudo')
subplot(2,3,5), imshow(img_liminvfilt), title('limited inv filt')


