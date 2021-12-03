% What is image enhancement?
% e.g. contrast enhance, deblur motion, eadge enhance, de-noises
% the goal: improve the "interpretability" of the images, or providing the
% "better inputs".

clear, clc
close all

% Original image
img = double(imread('tire.tif'));
figure, imagesc(img)
colormap(gray(256))
title('Original')

%% Power-law transform
c = 1;
gamma = .2;
s = c*img.^gamma;
figure, imagesc(s)
colormap(gray(256))
title('transformed')

%% Histogram equalization
% histogram: # of pixels with value in that range vs. bins
N = 64;
J = histeq(img, N); % N discrete level
figure, imagesc(J)
colormap(gray(256));
title('Histogram equalization')

%% Contrast limited adaptive histogram equalization (CLAHE)
J = adapthisteq(img, 'ClipLimit', 0.02);
% NumTiles: [row, col] to clip the image
% ClipLimit: higher number result in more contrast, default 0.01
figure, imagesc(J)
colormap(gray(256));
title('CLAHE') 

%% Laplacian
clear, clc
close all

img = double((imread('moon.tif')));

% Laplacian
kernal = [
    -1 -1 -1;...
    -1 8 -1;...
    -1 -1 -1];
L = conv2(img, kernal, 'same'); % 2nd derivative
L_disp = (L-min(L(:)))/(max(L(:))-min(L(:)))*255; % scaled for display purpose
img_L = img+L; % edge enhanced image

figure('Name', 'Laplacian'),
subplot(2,2,1), imshow(img, [0 255]), axis off
subplot(2,2,2), imshow(L, [0 255]), axis off
subplot(2,2,3), imshow(L_disp, [0 255]), axis off
subplot(2,2,4), imshow(img_L, [0 255]), axis off
colormap gray
