function J = region_growing(B)

B = logical(B);

% select the inital point to grow
X0 = false(size(B));
X0(170, 230) = true; 

% perform region filling
img_c = ~img; % complement of image
se = true(3); % structure element

% iter all
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

J = B;
J(Xk) = true;

end