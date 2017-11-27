function [depth] = getDepthFromNormals(n, mask)
  % [depth] = getDepthFromNormals(n, mask)
  %
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %
  imsize = size(mask);
  M = imsize(1);
  N = imsize(2);
  A = sparse(2*M*N, M*N);
  v = zeros(2*M*N, 1);
  nx = n(:,:,1);   % M x N matrix including x components of n
  ny = n(:,:,2);   % M x N matrix including y components of n
  nz = n(:,:,3);   % M x N matrix including z components of n
  
  % YOU NEED TO COMPLETE THIS.
  for j=1:N
      for i=1:M
        if mask(i,j) > 0
            idx = i + (j-1)*M;
            A( idx, idx   ) = -nz(i, j);
            A( idx, idx+M ) =  nz(i, j);
            A( idx + M*N, idx   ) = -nz(i, j);
            A( idx + M*N, idx+1 ) =  nz(i, j);
            v( idx )       = -nx(i, j);
            v( idx + M*N ) = -ny(i, j);
        end
      end
  end
  Z = A \ v;
  depth = reshape(Z, M, N);
  depth(mask==0) = 0;
  
  return