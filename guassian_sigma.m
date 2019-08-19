clc;
clear all;
close all;

% This script is used for find the optimal sigma for Gaussian filter.

%% load the original figure

fig = double(rgb2gray(imread('vesselTestImage.png')))/255;
[m,n] = size(fig);

%%

% Determine widths and standard deviations
width1 = 2; sigma1 = width1 / 2;
width2 = 10; sigma2 = width2 / 2;

test1 = imgaussfilt(fig,sigma1);
test2 = imgaussfilt(fig,sigma2);

figure();
subplot(3,1,1)
imagesc(fig);
title('a) original figure')
subplot(3,1,2)
imagesc(test1)
title(['b) Blurred image with \sigma = ' num2str(sigma1)]);
subplot(3,1,3)
imagesc(test2)
title(['c) Blurred image with \sigma = ' num2str(sigma2)]);

%% use quiver function to visualize the velocity vector
H1 = Hessian(test1);
for i = 1:m*n
[~,~,V1] = svd(H1(:,:,i));
v1(i) = V1(1,2);
v2(i) = V1(2,2);
end

v1 = reshape(v1,m,n)';
v2 = reshape(v2,m,n)';
[x,y] = meshgrid(1:1:m);

figure();
subplot(2,1,1)
imagesc(test2)
hold on
quiver(x,y,v1,v2)
title('a) velocity plot for \sigma = 1')

H2 = Hessian(test2);
for i = 1:m*n
[~,~,V2] = svd(H2(:,:,i));
v1(i) = V2(1,2);
v2(i) = V2(2,2);
end

v1 = reshape(v1,m,n)';
v2 = reshape(v2,m,n)';

subplot(2,1,2)
imagesc(test1)
hold on
quiver(x,y,v1,v2)
title('b) velocity plot for \sigma = 5')