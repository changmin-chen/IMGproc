%% Init
clc
close all
I = imread('afmsurf.tif');

%% parameters setting
% structuring element
% should design the se that has similar and larger shape compared to the target objects.
se = strel('disk', 15);

% H value for H-minima transform (H is a nonnegative scalar)
% minima transform detects all the intensity valleys "deeper" than a particular threshold H
H = 22;

%% Step 1: image enhancement
% enhanced image = (original + tophat) - bothat
I_en = imenhance(I, se);
figure,
subplot(1,2,1), imshow(I), title('original');
subplot(1,2,2), imshow(I_en), title('enhanced');

%% Step 2: reduce local minima (impose deep minima) of the image
Ic = imcomplement(I_en);
[I_impose, marker] = imimpose(Ic, H);
tmp  = I_en;
tmp(marker) = 255;
figure,
subplot(1,3,1), imshow(I_en), title('before minima-imposing');
subplot(1,3,2), imshow(tmp), title('overlayed marker image');
subplot(1,3,3), imshow(I_impose), title('after minima-imposing')

%% Step 3: watershed segmentation
wat = watershed(I_impose);
rgb = label2rgb(wat);
figure, 
imshow(rgb), title('watershed-segmented image')

%% helper functions
function I_en = imenhance(I, se)
% imtophat = I - imopen(I, se)
% returns image contains objects that are:
% "smaller" than the se and are "brighter" than their surrounding
% purpose: remove nonuniform-lighting background & remove objects larger than se.
% ---
% imbothat = imclose(I, se) - I
% returns image contains objects that are:
% "smaller" than the structuring element and are "darker" than their surrounding
% purpose: returns small spaces that are surrounding the objects
% ---
% image enhancement:
% I_en = (I + Itop) - Ibot
% ---
% se: structuring element, should design the se that has similar and larger shape
% compared to the target objects.
% imadd(A,B): return A+B, but would truncate exceeded values based on datatype
% imsubstract(A,B): return A-B, but negative result are rounded to 0.
Itop = imtophat(I, se);
Ibot = imbothat(I, se);
I_en = imsubtract(imadd(Itop, I), Ibot);

end

function [I_impose, marker] = imimpose(I, H)
% H-minima transform
% H is a nonnegative scalar
% detects all the intensity valleys deeper than a particular threshold H
% ---
% marker:
% the positions of the deep minima detected
% ---
% imimposemin: 
% removes all local minima in the original image except the marker positions
% deep minima values in marker positions are set to 0.
marker = imextendedmin(I, H);
I_impose = imimposemin(I, marker);

end