%% Section 1:  Basic image denoising
% using the following filters:
% (1) Arithmetic mean filter
% (2) Geometric mean filter
% (3) Harmoic mean filter
% (4) Contraharmonic mean filter
clear, clc
close all
img = double(imresize(rgb2gray(imread('saturn.png')), [512, 512]));
img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
img_n = imnoise(img, 'gaussian', 0, 0.05);

bk_sz = [5, 5];
% Arithmetic mean filter
% f(x,y) = mean(S,'all') for S in Processing Blocks
img_am = colfilt(img_n, bk_sz, 'sliding', @mean);

% Geometric mean filter
% f(x,y) = geomean(S,'all') for S in Processing Blocks
img_gm = colfilt(img_n, bk_sz, 'sliding', @geomean);

% Harmonic mean filter
% f(x,y) = harmmean(S, 'all') for S in Processing Blocks
img_hm = colfilt(img_n, bk_sz, 'sliding', @harmmean);

% Contraharmonic mean filter
% f(x,y) = charmmean(S, 'all') for S in Processing Blocks
img_chm = colfilt(img_n, bk_sz, 'sliding', @charmmean);

figure,
subplot(2,3,1), imshow(img), title('Original image')
subplot(2,3,2), imshow(img_n), title('Noisy image')
subplot(2,3,3), imshow(img_am), title('Arithmetic mean filter')
subplot(2,3,4), imshow(img_gm), title('Geometric mean filter')
subplot(2,3,5), imshow(img_hm), title('Harmonic mean filter')
subplot(2,3,6), imshow(img_chm), title('Contrarmonic mean filter')

%% Section 2: Adaptive denoising
% use MATLAB function "wiener2" to perform 
% 2D adaptive noise-removal filtering
clear, clc
close all
img = double(imresize(rgb2gray(imread('saturn.png')), [512, 512]));
img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
img_n = imnoise(img, 'gaussian', 0, 0.05);

% apply wiener2 function
bk_sz = [5, 5];
img_wie = wiener2(img_n, bk_sz);

figure,
subplot(1,3,1), imshow(img), title('Original image')
subplot(1,3,2), imshow(img_n), title('Noisy image')
subplot(1,3,3), imshow(img_wie), title('Adaptive mean filter')

%% Section 3: Adaptive median filter
% denoising Salt & Pepperer noise using adaptive median filter
clear, clc
close all
img = double(imresize(rgb2gray(imread('saturn.png')), [64, 64]));
img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
img_n = imnoise(img, 'salt & pepper', 0.05);
img_amf = adaptiveMedFilt(img_n);


figure,
subplot(1,3,1), imshow(img), title('Original image')
subplot(1,3,2), imshow(img_n), title('Noisy image')
subplot(1,3,3), imshow(img_amf), title('Adaptive median filter')

%% helper functions
function output = charmmean(C)
% calculate contra-harmonic mean
% for each block-column
% defalut order Q=2
Q = 2;
numerator = sum(C.^(Q+1), 1);
denominator = sum(C.^Q, 1);
output = numerator./denominator;

end

function J = adaptiveMedFilt(img)
% perform filtering using adaptive median filter
% input: image to filt
% output: filtered image
% latest modified date: 2021-12-08. Chen, Chang-Min

J = zeros(size(img));
for i = 1:size(img,1)
    for j = 1:size(img,2)
        bk_sz = [3, 3]; % initial block size
        
        % step 1: get inital block
        B = get_block(img, i, j, bk_sz);
        
        % step 2: level A evaluation
        go2B = levelA_eval(B);
        
        % step 3: increase the block size
        % until it pass the level A, and go to level B.
        % or until the block size equal the image size, and output the Zxy.
        while ~go2B
            bk_sz = bk_sz+2;
            B = get_block(img, i, j, bk_sz);
            go2B = levelA_eval(B);
            
            % if can't pass the level A, output the Zxy and break
            if size(B) == size(img)
                go2B = False;
                output = img(i, j);
                break
            end
        end
        
        % step 4: level B evaluation (should have pass the level A)
        % output the Zxy if pass the level B
        % output the median of the block if not
        if go2B  % check whether pass the level A
            B1 = img(i,j) - min(B(:));
            B2 = img(i,j) - max(B(:));
            if (B1>0) && (B2<0)
                output = img(i,j);
            else
                output = median(B(:));
            end
        end
        
        % step 5: assign the output
        J(i, j) = output;
        
    end
end


% local helper functions for Adaptive Median FIiter
%------------helper func 1-----------------
function B = get_block(img, i, j, bk_sz)
% get valid block with size bk_sz from img
% block center is [i, j]
ii = [i-(bk_sz(1)-1)/2, i+(bk_sz(1)-1)/2]; % [min i idx, max i idx]
jj = [j-(bk_sz(2)-1)/2, j+(bk_sz(2)-1)/2]; % [min j idx, max j idx]

if ii(1)<1 % validity min row
    ii(1) = 1;
end
if ii(2)>size(img,1) % validity max row
    ii(2) = size(img,1);
end
if jj(1)<1 % validity min column
    jj(1) = 1;
end
if  jj(2)>size(img,2) % validity max column
    jj(2) = size(img,2);
end
B = img(ii(1):ii(2), jj(1):jj(2));
end

%--------helper func 2---------
function chk = levelA_eval(B)
A1 = median(B(:)) - min(B(:));
A2 = median(B(:)) - max(B(:));
chk = (A1>0) & (A2<0);
end

end

