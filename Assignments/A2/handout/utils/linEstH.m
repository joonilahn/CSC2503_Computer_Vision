function [H, Sa] = linEstH(left, right, NUM_RESCALE)
  % [H, Sa, Sf] = linEstH(left, right, NUM_RESCALE)
  % Estimate the homography matrix H from two 3 x n matrices of
  % corresponding points left and right.  
  % Here left(:,k) x (H * right(:,k)) apprx= 0 where x denotes
  % the 3D cross product.
  % NUM_RESCALE (default TRUE) uses Hartley's rescaling. Always use
  % rescaling, unless you wish to show how badly the un-normalized
  % algorithm works.
  % Returns H along with the singular values Sa of the 2nPts x 9 homogeneous
  % linear system for H.
    
  if nargin < 3
    NUM_RESCALE = 1;
  end
  
  nPts = size(left,2);
  if nPts < 4 | nPts ~= size(right,2)
    fprintf(2,'lineEstH: Innappropriate number of left and right points.');
    H = [];
    return;
  end
  
  if size(left,1) == 2
    left = [left; ones(1, nPts)];
  else % Normalize to pixel coords
    left = left./repmat(left(3,:), 3,1);
  end
  if size(right,1) == 2
    right = [right; ones(1, nPts)];
  else % Normalize to pixel coords
    right = right./repmat(right(3,:), 3,1);
  end
  
  imPts = cat(3, left, right);
  
  %% Rescale image data for numerical stability.
  if NUM_RESCALE
    Knum = repmat(eye(3), [1,1,2]);
    %%% Rescale for numerical stability
    mn = sum(imPts(1:2,:,:),2)/nPts;
    mns = reshape(mn, [2 1 2]);
    var = sum(sum((imPts(1:2,:,:)-repmat(mns, [1 nPts 1])).^2,2)/nPts, 1);
    %% Scale image points so that sum of variances of x and y = 2.
    scl = sqrt(2./var(:));
    %% Sanity: varScl =  var .* reshape(scl.^2, [1 1 2]) % Should be 2
    
    %% Scale so x and y variance is roughly 1, translate so image mean (x,y) is zero.
    Knum(1:2,3,:) = -mn;
    Knum(1:2,:,:) = Knum(1:2,:,:).*repmat(reshape(scl, [1 1 2]), [2, 3,1]);
    for kIm = 1:2
      imPts(:,:,kIm) = reshape(Knum(:,:,kIm),3,3) * imPts(:,:,kIm);
    end
    %% Sanity check
    % sum(imPts(1:2,:,:),2)/nPts  % Should be [0 0]'
    % sum(sum(imPts(1:2,:,:).^2,2)/nPts,1) % Should be 2.
  end

  %% Make constraint matrix A. 
  %% You fill in below here.
  left = reshape(imPts(:,:,1), [3 nPts]);   % 3xK matrix
  right = reshape(imPts(:,:,2), [3 nPts]);  % 3xK matrix
  A = [];          % left_x: (3*nPts) x 3 cross product matrix
  for k=1:nPts
    left_k = left(:,k);
    right_k = right(:,k);
    A_k = [0 0 0 -1 -1 -1   left_k(2)  left_k(2)  left_k(2);
           1 1 1  0  0  0  -left_k(1) -left_k(1) -left_k(1)];
    A_k = A_k .* repmat(right_k',2,3);
    A = [A; A_k];
  end
  
  %% Factor A
  [Ua Sa Va] = svd(A); Sa = diag(Sa);
    
  %% Set H to be the right null vector of A, reshaped to a 3x3 matrix.
  H0 = reshape(Va(:,end), 3,3)';
  
  %% Undo the renormalization
  if NUM_RESCALE
    H = inv(reshape(Knum(:,:,1),3,3)) * H0 * reshape(Knum(:,:,2),3,3);
  end
  H = H ./ sqrt(sum(sum(H.^2)));
  
%   %% Sanity check
%   obj_val = 0;
%   for k=1:nPts
%       A_k = reshape(A(k,:),3,3);
%       check_sum = A_k * H * right(:,k);
%       obj_val = obj_val + sum(check_sum(1:2).^2);
%   end
%   fprintf('Objective value is %f\n', obj_val);
end  
  
  