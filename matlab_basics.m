%% Programming for performance
% Preallocation 先矩陣宣告大小，不要持續去改變矩陣大小
x = zeros(1, 100);
for i = 1: 100
    x(i+1) = x(i)+3;
end

% Vectorization 向量化
t = 0: .1: 10;
a = sin(t);

% Avoid global variables
% Avoid overloading build-in func
% Avoid using "data as code"

% Generating C code from MATLAB code
% Using MATLAB Workers on multi-core processors and clusters
% Using GPU

% Assign: row-major
for j = 1:n % cols
    for i = 1:n % rows
        x(i,j) = 3;
    end
end

% Avoid using if in loop
for i = 1: 100
    if i == 1
        % do things when i=1
    else
        % do things when i=2:100
    end
end

% Use func. with "real" form
x = 4;
x = sqrt(x);
x = realsqrt(x);

% Summation: row-major summation
A = rand(1000);
tic
for t=1:1000
    sum(A,2);
end
toc

tic
B = A';
for t=1:1000
    sum(B,1);
end
toc

%% use GPU
A = gpuArray(a); % Put a from RAM to GPU
B = gpuArray(b); % Put b from RAM to GPU
C = A*B;
c = gather(C); % Catch data from GPU

%% 
