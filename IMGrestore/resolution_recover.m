clear, clc
close all
% parameters setting
motion_length = 9;
theta = 45;

% load image
img = double(imread('cameraman.tif'));
img = (img-min(img(:)))/(max(img(:))-min(img(:)));

% add motion
kernel = fspecial('motion', motion_length, theta);
% img_m = imfilter(img, kernel); % motion

% add Gaussian
% img_m_b = imnoise(img_m, 'gaussian', 0, 0.02);

% display
% figure,
% subplot(1,3,1), imshow(img)
% subplot(1,3,2), imshow(img_m)
% subplot(1,3,3), imshow(img_m_b)