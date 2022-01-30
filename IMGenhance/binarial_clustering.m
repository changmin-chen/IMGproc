function clusters = binarial_clustering(BW)
% -----
% Example:
% clusters = binarial_clustering(BW);
% figure,
% for k = 1:size(clusters,3)
% imshow(clusters(:,:,k))
% title(['Cluster #' num2str(k)])
% pause(0.5)
% end
% -----

BW = logical(BW);
ndim = ndims(BW);
clusters = zeros(size(BW));
ncluster = 0;
tic;
while sum(sum(BW))>0
    ncluster = ncluster +1;
    new_cluster = region_growing(BW);
    if ndim==2
        clusters(:,:,ncluster) = new_cluster;
    elseif ndim==3
        clusters(:,:,:,ncluster) = new_cluster;
    end  
    BW(new_cluster) = false;
    
end
fprintf('Clustering completed.\n')
toc;
fprintf([num2str(ncluster), ' clusters were found.\n'])

end

function Xk = region_growing(BW)
% -----
% region_growing(BW)
% -----
% input:
% BW: binary image, could be 2 or 3 dimensional
% -----
% neighborhood type: face + edge
% -----

% 2D neighborhood definition
neighb_2d = true(3);

% 3D neighborhood definition
block = false(3);
block(2,2) = true;
neighb_3d = cat(3, block, neighb_2d, block);

% select neighborhood base on input BW dimension
ndim = ndims(BW);
if ndim == 2
    neighb = neighb_2d;
elseif  ndim == 3
    neighb = neighb_3d;
else
    error('Input should be 2 or 3 dimensional.');
end

% select seed point, and define Xk & Xk+1
Xk = false(size(BW));
Xk(find(BW>0, 1, 'first')) = true;
Xkk = imdilate(Xk, neighb) & BW;

% loop
iteration_count = 1;
while ~all(Xkk==Xk, 'all') % if not all Xk+1 equal to Xk
    Xk = Xkk;
    Xkk = imdilate(Xk, neighb) & BW;
    iteration_count = iteration_count + 1;
end

end
