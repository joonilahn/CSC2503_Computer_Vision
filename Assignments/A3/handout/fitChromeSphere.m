function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
    
  mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
  mask = mask(:,:,1) / 255.0;

  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,n) = im(:,:,1);           % red channel
  end

  % YOU NEED TO COMPLETE THIS FUNCTION

  % Compute center and radius of the sphere
  [sphere_y sphere_x] = find(mask>0);    % find pixels representing sphere
  center_x = (max(sphere_x) + min(sphere_x)) / 2;  % center_x
  center_y = (max(sphere_y) + min(sphere_y)) / 2;  % center_y
  center = [center_x; center_y; 0];                % center_z = 0
  radius = ((max(sphere_x) - center_x) + (max(sphere_y) - center_y)) / 2;
  
  % loop over the images to estimate the direction of the light source
  d_c = [0; 0; -1];     % the direction of camera
  L = zeros(3, nDir);   % 3 x nDir image of light source directions
  for i=1:nDir
      % Find the brighestest point
      [argmax_y argmax_x] = find(imData(:,:,i)==max(max(imData(:,:,i))));
      X = [median(argmax_x); median(argmax_y)] - center(1:2);
                            % The brightest point X
      X_z = -sqrt(radius^2 - X(1)^2 - X(2)^2);
      X = [X; X_z];         % Append z-component
      normal = X / norm(X);  
                            % Normal vector at X
      L(:,i) = 2 * dot(normal, d_c) * normal - d_c;  % Direction of reflection
  end
  
  return;

