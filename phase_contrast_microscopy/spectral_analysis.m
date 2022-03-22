%% init
% img1: original image
% img2: placed a filter at Fourier plane, which was a "low-pass filtered image".
clear, clc
close all
img1 = rgb2gray(imread('IMG_3965.JPG'));
img2 = rgb2gray(imread('IMG_3966.JPG'));

%% spatial domain analysis
figure, imshow(img1, []), title('oringal image')
figure, imshow(img2, []), title('low-pass filtered image'),
figure, imshow(img1-img2, []), title('difference b/w two image')
figure,
histogram(img1, 32)
hold on
histogram(img2, 32)
legend({'original', 'low-pass filtered'})
xlabel('grayscale pixel values')
xlim([0 255]) 
xticks(0:25:255)
ylabel('number of pixels')
title('Image histogram')

% compute energy
energy1 = sum(sum(img1.^2));
energy2 = sum(sum(img2.^2));
fprintf(['Energy for image1= ' num2str(energy1), '.\n'...
    ,'Energy for image2= ' num2str(energy2) '.\n'])
% nearly same level of energy

% compute noise level
roi1 = img1(50:450, 1000:1500);
roi2 = img2(50:450, 1000:1500);
figure,
subplot(2,2,1), imshow(roi1, []), title('ROI in original image')
subplot(2,2,2), histogram(roi1, 64), title('histogram original')
subplot(2,2,3), imshow(roi2, []), title('ROI in low-pass filtered image')
subplot(2,2,4), histogram(roi2, 64), title('histogram low-pass filtered')

noise1 = std(double(roi1(:)));
noise2 = std(double(roi2(:)));
fprintf(['Noise level for image1= ' num2str(noise1), '.\n'...
    ,'Noise level for image2= ' num2str(noise2) '.\n'])
% nearly same level of noise

%% frequency domain analysis
f1 = fft2(img1);
f1shift = fftshift(f1);
mag_spectrum1 = 10*log10(abs(f1shift));

f2 = fft2(img2);
f2shift = fftshift(f2);
mag_spectrum2 = 10*log10(abs(f2shift));

figure,
imshow(mag_spectrum1, []), title('original spectrum (dB)')
colormap jet
colorbar

figure,
imshow(mag_spectrum2, []), title('low-pass filtered spectrum (dB)')
colormap jet
colorbar

%% continue. Comparison in frequency domain
mag_spectrum_diff = mag_spectrum1-mag_spectrum2;
mag_spectrum_diff_BW = mag_spectrum_diff>0;
figure,
imshow(mag_spectrum_diff_BW), title('original spectrum > low-pass filtered spectrum')
colormap([0 0 1; 1 1 0])
colorbar
% conclusion:
% low-pass filter was not ideally placed at the Fourier plane

%% simulate filtering
% draw ellipse
sz = size(img1);
a = 1000;
b = 3000;
ep1 = draw_ellipse(sz, 1000, 3000);
ep2 = draw_ellipse(sz, 3000, 1000);
ep = ep1| ep2;
figure,
imshow(~ep)

% draw circle
cir = draw_circle(sz, 150);
figure,
imshow(cir)

% extract high-freq component from img1
f1shift_high = f1shift;
f1shift_high(ep) = 0;
img1_high = ifft2(ifftshift(f1shift_high));
img1_high = abs(img1_high);
figure, imshow(img1_high, [])

% extract low-freq component from img1
f1shift_low = f1shift;
f1shift_low(~cir) = 0;
img1_low = ifft2(ifftshift(f1shift_low));
img1_low = abs(img1_low);
figure, imshow(img1_low, [])

function circle = draw_circle(sz, radius)
circle = false(sz);
center = round(sz/2);
[xx, yy] = meshgrid(1:sz(2), 1:sz(1));
d2 = (yy-center(1)).^2+(xx-center(2)) .^2;
circle(d2<radius^2) = true;

end

function ellipse = draw_ellipse(sz, a, b)
% ellipse = false(sz);
center = round(sz/2);
[xx, yy] = meshgrid(1:sz(2), 1:sz(1));
ellipse = (yy-center(1)).^2 ./ a^2 + (xx-center(2)).^2 ./b^2 <= 1;

end
