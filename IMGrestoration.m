clear, clc
close all

img = double(imread('eight.tif'));
img = (img-min(img(:))) / (max(img(:))-min(img(:)));

% Type 1: Gaussian noise
J = imnoise(img, 'gaussian', 0, 0.01);

% Type 2: Poisson noise
% img = img/1e12;
% J = imnoise(img, 'poisson');
% J = J*1e12;
% img = img*1e12;

% probe
roi_ori = img(1:80, 1:80);
h_ori = imhist(roi_ori, 256);
roi_noise = J(1:80, 1:80);
h_noise = imhist(roi_noise, 256);

figure,
subplot(3,2,1), imshow(img), title('original image')
subplot(3,2,2), imshow(J), title('noise image')
subplot(3,2,3), imshow(roi_ori), title('original background roi')
subplot(3,2,4), plot(linspace(0,1,256),h_ori), title('original histogram')
subplot(3,2,5), imshow(roi_noise), title('noised background roi')
subplot(3,2,6), plot(linspace(0,1,256),h_noise), title('noised histogram')

%% De-noise using wiener filter
clear, clc
close all
RGB = imread('saturn.png');
img = double(rgb2gray(RGB));
img = (img-min(img(:))) / (max(img(:))-min(img(:)));
J = imnoise(img, 'gaussian', 0, 0.025);
K = wiener2(J, [5,5]);

figure,
subplot(1,3,1), imshow(img), title('original image')
subplot(1,3,2), imshow(J), title('noised image')
subplot(1,3,3), imshow(K), title('de-noised image')

%% Adaptive median filter
% noise type: salt & papper noise


%%






