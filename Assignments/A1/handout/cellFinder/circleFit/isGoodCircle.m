function [goodCircle, klDiv] = isGoodCircle(x0, r, w, ...
                                     circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.
  
  x0 = x0(:);  % Decide, row or column
  
  % YOU FINISH THIS
  wSum = sum(w);                        % sum of all weights
  wMax = max(w);                        % max weight
  wThreshold = 0.5;                     % Threshold for weight: set 0.5
  strongW= w(find(w > wThreshold));     % Strong weights above threshold
  supportRatio = sum(strongW) / wSum;   % Ratio: sum(strongweights)/sum(allW)
  ratioThreshold = 0.7;                 % The ratio above may exceed this threshold
  
  if supportRatio > ratioThreshold
      goodCircle = true;
  else
      goodCircle = false;
  end
  
%   figure(5); bar(w);
end