function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.
  
x0 = initx0(:);  % Make it a column vector.

% FINISH THIS CODE
r = initr;                      % initial radius
K = size(pts, 1);               % number of points
a = -2 * x0;                    % a = -2*Xc
b = x0' * x0 - r^2;             % b = x_0'x_0 - r^2
e = zeros(K, 1);                % errors: Kx1 array
w = zeros(K, 1);                % weights: Kx1 array
psi = zeros(K, 1);              % psi(influence func): Kx1 array
prev_a = inf* ones(size(a));    % previous a: initially set to inf
prev_b = inf;                   % previous b: initially set to inf

% repeat until the growth a, b are significantly slow
while (sum((a-prev_a).^2) > 1e-3) && (sum((b-prev_b).^2) > 1e-3)
    e = pts * a + b + dot(pts, pts, 2); % error: a'x_i + b + x_i'x_i
    w = (2 * (sigmaGM^2)) ./ ((sigmaGM^2 + e.^2).^2); % weight=(1/e)*(drho/de)
    C = [pts ones(K,1)];                % C: Kx3 matrix
    d = sum(pts.^2,2);                  % d: Kx1 vector
    A = C' * diag(w) * C;                  % A: C'*diag(w)*C
    B = - C' * diag(w) * d;                % B: -C'*diag(w)*d
    sol = pinv(A) * B;                  % p = pseudoinv(A)*B
    prev_a = a;                         % a = prev_a
    prev_b = b;                         % b = prev_b
    a = sol(1:2);                       % a = p(1:2)
    b = sol(3);                         % b = p(3)
end

x0 = - a / 2;                           % xc = -a / 2
r = sqrt(x0' * x0 - b);                 % r = sqrt(x0'*x0 - b)
maxW = max(w);

end