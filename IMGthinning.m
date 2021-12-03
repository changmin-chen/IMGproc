clear, clc
close all

threshold = 160;
img = imread('ori_1.jpg');
img = 1 - (img>threshold);

figure,
subplot(1,2,1), imshow(img), title('original image')

%%
close all
img = rand(200, 300) > 0.3;

figure,
subplot(1,2,1), imshow(img), title('ori')

se1 = true(1,2);
se2 = true(1,3);
se3 = true(1,4);
se4 = true(1,5);
se5 = [0 0 0; 1 1 0; 1 1 0];
se6 = [1 1 1; 1 1 1; 1 1 1];
img_thin = imthin(img, se1, se2, se3, se4, se5, se6);
subplot(1,2,2), imshow(img_thin), title('thin image')

%%
a = hit_or_miss(img,se6);
imagesc(a)


%% helper functions

function SEout = padc_se(SEin)
% get padded complementary se
% SEin = X, SEout = W-X, and size(SEout) = size(SEin)+2
% pad two rows and two columns at first and last position
SEout = true(size(SEin)+2);
SEout(2:end-1, 2:end-1) = ~SEin;
end

function Bout = hit_or_miss(Bin, se)
% hit-or-miss processing for binary image Bin
% only the pixels that exactly matched the se preserved
Bin = logical(Bin);
se = logical(se);
A1 = imerode(Bin, se);
A2 = imerode(~Bin, padc_se(se));
Bout = A1 & A2;
end

function B = imthin(varargin)
B = logical(varargin{1});
ses = varargin(2:end);

for iter = 1: 10
    Bin_count = sum(sum(B));
    disp(num2str(Bin_count));
    for i = 1: length(ses)
        se = logical(ses{i});
        for j = 1:4
%             disp(num2str(sum(sum(B))));
            B_hm = hit_or_miss(B, se);
            B = logical(B - B_hm);
            se = rot90(se);
        end
        Bout_count = sum(sum(B));
    end
    disp(num2str(Bout_count));
    if Bin_count == Bout_count, break, end
end

end