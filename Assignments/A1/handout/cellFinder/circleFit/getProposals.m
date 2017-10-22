function [circles] = getProposals(normals, p, numGuesses)
  % [circles] = getProposals(normals, p, numGuesses)
  % Attempt to produce up to numGuesses circle proposals from
  % the edgel data p and normals.  For typical data sets
  % we will be able to produce numGuesses proposals.  However,
  % on some datasets, say with only a few edgels, we may not be
  % able to generate any proposals.  In this case size(circles,1)
  % can be zero.
  % Input:
  %  normals - N x 2 edgel normals
  %  p         N x 2 edgel positions
  %  numGuesses - attempt to propose this number of circles.
  % Return:
  %   circles a P x 3 array, each row contains [x0(1) x0(2) r]
  %           with 0 <= P <= numGuesses.
  
  
  % YOU NEED TO FILL IN CODE HERE.
  if numGuesses < size(p,1)
      P = numGuesses;
  else
      P = size(p,1);
  end
  
  circles = zeros(P,3);     % generate initial P x 3 array
  N = size(p,1);            % number of edgels
  sampleIdx = randperm(P);  % Px1 array with random numbers
  for i=sampleIdx
      x1 = p(i, :);        % x1: first random sample
      n1 = normals(i, :);  % n1: x1's normal vector
      n1DotNormals = normals * n1'; % dot product (n1, All normals)
      minAngleIdx = find(abs(n1DotNormals) < 0.01);  % Find perpendicular
                                               % (90 degree) normal vector
      
      if size(minAngleIdx,1) == 0           % if not found,
        [m minAngleIdx] = min(abs(n1DotNormals));   % find the closest 
                                                    % angle to 90 degree
      end
      if size(minAngleIdx,1) > 1            % if found more than one
        candidates = p(minAngleIdx, :);
        dist = [];
        for j=1:size(candidates,1)
            dist(j) = sum((x1 - candidates(j)).^2);
            %dist = sqrt(sum((x1 - candidates).^2, 2)); % pick the closest point
        end
        [m minAngleIdx] = min(dist);
      end
      
      x2 = p(minAngleIdx, :);         % x2's normal vector might be
                                      % perpendicular to the x1's
      n2 = normals(minAngleIdx, :);
      theta = acos(dot(n1, n2));      % calculate the angle between n1, n2
      dist = sqrt(sum((x1-x2).^2));   % dist: distance between x1, x2
      radius = 0.5 * dist * abs(cos(pi-theta/2));
                                % radius = 0.5 * distance * cos(pi-theta/2)
      center = x1 + radius * n1;      % center = x1 + radius * n1(direction)
      circles(i,:) = [center radius];
  end
  
%   Visualize
%   for i=1:size(circles,1)
%     figure(1); hold on;
%     plotCircle(circles(i, 1:2)', circles(i,3), 'r');
%   end
%   
end