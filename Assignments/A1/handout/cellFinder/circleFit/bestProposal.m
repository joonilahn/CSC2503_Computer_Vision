function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals
  
  % YOU COMPLETE THIS
  K = size(circles,1);      % number of circles
  P = size(p,1);            % number of edgels
  errors = zeros(K,1);      % errors for circles
  for i=1:K
    xc_i = circles(i,1:2);  % ith circle
                            % compute error_i = ||p_i-xc_i||^2 - r^2
    error = [];
    for j=1:P
        error(j) = norm(p(j,:) - xc_i)^2 - circles(i,3)^2;
    end
    negError = find(error<0);
    error(negError) = error(negError).^2;
    error = sum(error);     % sum of error for ith circle
    errors(i) = error;
  end
  
  [min_error min_idx] = min(errors);    % find minimum error and its index  
  circle = circles(min_idx, :);         % select the circle with min error
  
  % Visulize
%   figure(1); hold on;
%   plotCircle(circle(1:2)', circle(3), 'b');

end
  