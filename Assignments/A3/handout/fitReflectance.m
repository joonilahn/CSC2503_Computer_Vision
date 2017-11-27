function [n, albedo] = fitReflectance(im, L)
  % [n, albedo] = fitReflectance(im, L)
  % 
  % Input:
  %   im - nPix x nDirChrome array of brightnesses,
  %   L  - 3 x nDirChrome array of light source directions.
  % Output:
  %   n - nPix x 3 array of surface normals, with n(k,1:3) = (nx, ny, nz)
  %       at the k-th pixel.
  %   albedo - nPix x 1 array of estimated albdedos
    

  % YOU NEED TO COMPLETE THIS
  g = im * L' * inv(L*L');      % solve GLL' = IL'
  albedo = sqrt(sum(g.^2, 2));  % magnitude of g
  n = g ./ repmat(albedo, 1, 3) % direction of g
  
  % Find NaN and set to be zero
  % There are NaN values in n if there are zero input in im matrix.
  % We cannot estimate n vector for the pixels
  nanidx = find(isnan(n));
  n(nanidx) = 0;
  
  return;


