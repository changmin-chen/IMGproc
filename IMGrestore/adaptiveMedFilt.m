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

end

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

function chk = levelA_eval(B)
A1 = median(B(:)) - min(B(:));
A2 = median(B(:)) - max(B(:));
chk = (A1>0) & (A2<0);

end
