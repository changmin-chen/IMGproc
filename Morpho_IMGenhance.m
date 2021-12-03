%% Fingerprint
clear, clc
close all

threshold = 160;
img = imread('ori_1.jpg');
img = img<=threshold;

% se
img_proc = img;
se1 = strel('square', 2);
se2 = strel('line', 1, 135);

% erode dilate = open
img_proc = imerode(img_proc, se1);
img_proc = imdilate(img_proc, se1);

% dilate erode = close
img_proc = imdilate(img_proc, se2);
img_proc = imerode(img_proc, se2);

% equivalent to the above
% img_proc = imopen(img_proc, se1);
% img_proc = imclose(img_proc, se2);

figure, subplot(1,2,1), imagesc(img), axis off
subplot(1,2,2), imagesc(img_proc), axis off
colormap(gray)

%% Hit or miss
clear, clc
close all
img = imread('test.tiff');

%% Region filling
clear, clc
close all
load('test.mat');



