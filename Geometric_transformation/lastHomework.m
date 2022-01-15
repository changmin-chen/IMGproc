%% Section 1: Image Fusion: Weighted averaging pixel-wise
clear, clc
close all
transparency = 0.1;

% registration
load('img_reg_fus.mat');
[optimizer, metric] = imregconfig('monomodal');
img2_reg = imregister(img2, img1,...
    'rigid', optimizer, metric);

% fusion
w = 0.7;
n_cback = 256;
n_cfore = 64;
% scale pixel values to 0~1
img1 = mat2gray(img1, [0, max(img1(:))]);
img2_reg = mat2gray(img2_reg, [0, max(img2_reg(:))]);
% adjust pixel values to the desired color indexes
img1 = img1*(n_cback-1)+1; % 256 steps, for index 1~256
img2_reg = img2_reg*(n_cfore-1)+n_cback; % 64 steps, for index 256~319
cmap = [gray(n_cback); jet(n_cfore)];
% specify AlphaData for the foreground image
alpha = ones(size(img2_reg));
alpha(img2_reg==min(img2_reg(:))) = 0;
alpha = alpha*transparency;

figure,
subplot(1,2,1), imshow(img1,[]), title('background')
subplot(1,2,2), imshow(img2_reg,[]), title('foreground')
figure,
image(img1)
hold on
image(img2_reg, 'AlphaData', alpha)
title(['Fused, foreground transparency = ' num2str(transparency)])
colormap(cmap)
colorbar
axis off

%% Section 2: Image Fusion: Wavelet
% Objective: fuse the multifocus images
% Reference: "Multifocus image fusion in wavelet domain, 2003"

clear, clc
close all
wname = 'db5';

% discrete 2d wavelet transform
img1 = double(imread('book1.jpg'));
img2 = double(imread('book2.jpg'));
[A1, H1, V1, D1] = dwt2(img1, wname, 'mode', 'per');
[A2, H2, V2, D2] = dwt2(img2, wname, 'mode', 'per');

% take the average of the baseband
A = (A1+A2)./2;

% take max values of high-frequency components
H = max(H1, H2);
V = max(V1, V2);
D = max(D1, D2);
fusion = idwt2(A, H, V, D, wname);

% display
figure('Units', 'normalized', 'Position', [0.02, 0.3, 0.15, 0.15], 'Name', 'clear background')
imshow(img1, [], 'border', 'tight'), 
figure('Units', 'normalized', 'Position', [0.4, 0.3, 0.15, 0.15], 'Name', 'clear foreground')
imshow(img2, [], 'border', 'tight')
figure('Units', 'normalized', 'Position', [0.7, 0.3, 0.15, 0.15], 'Name', 'fused image')
imshow(fusion, [], 'border', 'tight')
