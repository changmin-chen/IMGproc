%% Section 1.1: Unsharp filtering
% formula: f_sharp = f + w*(f - f_low-pass)
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

%% Section 1.2: High boost filtering
% formula: f_sharp = (A-1)*f + f_high-pass
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

%% Section 1.3: High frequency emphasis filter
% formula: H_hfe = a + b*H_hp
% a>=0 and b>a
% a = 0.25~0.5 and b=1.5~2.0
% functions:
% (1) fspecial: to create the filter
% (2) imfilter: to filter the image

clear, clc
close all
img = imread('tire.tif');

a = 0.5;
b = 1.5;
filter = fspecial('laplacian'); % a high-pass filter
img_hp = imfilter(img, filter);
img_hfe = a + b*img_hp;
img_hfe = imhistmatch(img_hfe, img);

figure,
subplot(1,3,1), imshow(img), title('original image')
subplot(1,3,2), imshow(img_hp), title('high-pass filtered image')
subplot(1,3,3), imshow(img_hfe), title('high freq enphasis filtered')

%% Section 2: Process fingerprint image
% using basic morphological processing technique
% including: erosion, dilation, close, open
% dilation then erosion = close
% erosion then dilation = open

clear, clc
close all

threshold = 160;
img = imread('ori_1.jpg');
img = 1 - (img>threshold);

% morphological processing
se1 = [0 1 0; 0 1 0; 0 1 0];
img_d1 = imdilate(img, se1);
img_ed1 = imerode(img_d1, se1);
bg1 = ~img_ed1;

se2 = se1';
img_d2 = imdilate(img, se2);
img_ed2 = imerode(img_d2, se2);
bg2 = ~img_ed2;

img_proc = ~(bg1 & bg2);

figure,
subplot(2,3,1), imshow(img), title('original image')
subplot(2,3,2), imshow(img_d1), title('image dilated by se1')
subplot(2,3,3), imshow(img_ed1), title('image closed by se1')
subplot(2,3,4), imshow(img_d2), title('image dilated by se2')
subplot(2,3,5), imshow(img_ed2), title('image closed by se2')
subplot(2,3,6), imshow(img_proc), title('processed image')
% the processed image may be more connected compared to original image
% i.e. less holes then original. However, it losed details

%% Section 3: Hit-or-miss
% a basic morphological tool for shape detection
% the pixels lasted were the locations of the regions that exactly matched the shape X
clear, clc
close all
img = imread('test.tiff');

% define the target shape X
X = [0 1 0; 1 1 1; 0 1 0]; % X, target shape
W = ones(5); 
W(2:end-1, 2:end-1) = W(2:end-1, 2:end-1) - X; % W-X, the window

% peform hit-or-miss algorithm
img_er = imerode(img, X);
img_er_cm = imerode(~img, W);
img_hm = img_er & img_er_cm;

figure,
subplot(2,3,1), imshow(img), title('original image')
subplot(2,3,2), imshow(X), title('X, target shape (size 3x3)')
subplot(2,3,3), imshow(W), title('W-X, window')
subplot(2,3,4), imshow(img_er), title('image eroded by X')
subplot(2,3,5), imshow(img_er_cm), title('image^c eroded by W-X')
subplot(2,3,6), imshow(img_hm), title('image after hit-or-miss')

%% Section 4: Region filling
clear, clc
close all
load('test.mat');
img = logical(img);

% select the inital point to grow
X0 = false(size(img));
X0(170, 230) = true; 

% perform region filling
img_c = ~img; % complement of image
se = true(3); % structure element

% demo iter 1
X1 = imdilate(X0, se) & img_c;

% demo iter 10
X10 = X0;
for k = 1: 10
    X10 = imdilate(X10, se) & img_c;
end

% demo iter all
Xk = X0;
tic;
Xkk = imdilate(Xk, se) & img_c; % Xk+1
iter = 1;
while ~all(Xkk==Xk, 'all') % if not all Xk+1 equal to Xk...
    Xk = Xkk;
    Xkk = imdilate(Xk, se) & img_c;
    iter = iter+1;
end
disp(['Region filling completed after ' num2str(iter) ' iteration.']);
toc;

img_filled = img;
img_filled(Xk) = true;

figure,
subplot(2,3,1), imshow(img), title('original image')
subplot(2,3,2), imshow(X0), title('X0, initial pixel to grow')
subplot(2,3,3), imshow(X1), title('X1, dilate after 1 iteration')
subplot(2,3,4), imshow(X10), title('X10, dilate after 10 iteration')
subplot(2,3,5), imshow(Xkk), title('Xk, dilate until Xk=Xk-1')
subplot(2,3,6), imshow(img_filled), title('filled image')
