% File: dinoTrueF.m
% Assignment 2 help script

%% Get estimated F from dinoTestF
dinoTestF
F1 = F;

%% Compute true fundamental matrix
Tvec = dLeft - dRight;
T = zeros(3,3);
T(1)=0;
T(2) = Tvec(3);
T(3) = Tvec(2);
T(4) = Tvec(3);
T(5) = 0;
T(6) = -Tvec(1);
T(7)=-Tvec(2);
T(8)=Tvec(1);
T(9)=0;
E = Rleft* T * Rright';                         % Essential matrix
F0 = pinv(MintLeft)' * E * pinv(MintRight);     % Ground-truth F matrix


%% A.2. Compare F0 and F1
hFig = figure(3);
clf; hold on;
% Plot all interest point locations as blue .'s
[X, Y] = meshgrid(xbins, ybins);
plot(X, Y, '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image');

nPoints = 8;  %% Number of points per axis
nCol = 16;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.
xbins = linspace(ax(1), ax(2), nPoints);
ybins = linspace(ax(2), ax(4), nPoints);

for i = 1:nPoints
    for j = 1:nPoints
      % Plot interest point location corresponding to epipolar line as a "o"
      % in the same colour as the epipolar line.
      plot(xbins(i), ybins(j), 'o', 'Color', col(mod(k,nCol)+1,:));
      % Plot epipolar line.
      lk = imPts(:,k,2)' * F';
      epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
      set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
    end
end