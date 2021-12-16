%% Draw circle with for loop
clear, clc
close all
map = zeros(50); % background
c = round(size(map)/2); % center (x,y)
radius = 10;

for j = 1:size(map,1)
    for i = 1:size(map,2)
        if  sqrt(sum(([i, j]-c).^2)) >= radius % distance = sqrt((x-xc)^2+(y-yc)^2)
            map(i,j) = 0;
        end
    end
end
f = figure;
a = axes;
imshow(map,[], 'parent' ,a)
title('Draw circle with for loop')
axis(a,'image')

%% Draw circle without for loop
clear, clc
close all
map = zeros(50); % background
c = round(size(map)/2); % center (x,y)
radius = 10; % circle radius

% distance = sqrt((x-xc)^2+(y-yc)^2)
[x, y] = meshgrid(1:size(map,1), 1:size(map,2)); % all (x,y) pair
d = sqrt(((x-c(1)).^2+(y-c(2)).^2)); % distance to circle center
map(d<=radius) = 1; % circle = distance<=ridus

f = figure;
a = axes;
imshow(map,[], 'parent' ,a)
axis(a,'image')

%% Vectorization
clear, clc
n = 10; m = 10;
jb = 2;
je = n+1;
ib = 3;
ie = m+1;
u = rand(ie,je);
omega = 0.5;
del = 0;

% original code format
for i = ib:2:n
    for j = jb:2:m
        up = (u(i,j+1)+u(i+1,j)+u(i-1,j)+u(i,j-1))*0.25;
        u(i,j) = (1-omega)*u(i,j)+omega*up;
        del = del + abs(up-u(i,j))
    end
end

% vectorized code format
i = ib:2:n;
j = jb:2:m;
up = (u(i,j+1)+u(i+1,j)+u(i-1,j)+u(i,j-1))*0.25;
u = (1-omega).*u + omega.*up;
% conv
del2 = sum(sum(abs(up-u)))
