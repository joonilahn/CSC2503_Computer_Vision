% File: dinoTrueF.m
% Assignment 2 help script

%% Get estimated F from dinoTestF
dinoTestF
close(1); close(2);
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
nPoints = 8;  %% Number of points per axis
nCol = 16;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

figure(1);
xbins = linspace(0, 100, nPoints);
ybins = linspace(0, 100, nPoints);
[X, Y] = meshgrid(xbins, ybins);
plot(X, Y, '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
% axis([0 100 0 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Regularly spaced grid of points');

hfig = figure(2);
clf; hold on;
title('Epipolar Lines for the grid of points');
axis([-150 150 -100 100]);
% axis([0 100 0 100]); 
% axis xy; axis equal;
perpErr = [];
for i = 1:nPoints
    for j = 1:nPoints
      % Plot interest point location corresponding to epipolar line as a "o"
      % in the same colour as the epipolar line.
      plot(xbins(i), ybins(j), 'o', 'Color', 'b');
      % Plot epipolar line.
      lF0 = [xbins(i); ybins(j); 1]' * F0;
      lF1 = [xbins(i); ybins(j); 1]' * F1;
      ep0 = cropLineInBox(lF0(1:2), lF0(3), cropBox);
      if sum(sum(isnan(ep0))) > 0
          break;
      else
        set(line(ep0(:,1), ep0(:,2)), 'Color', 'r');
        maxperpErr_ij = max(abs([ep0 ones(2,1)]*lF1') / norm(lF1(1:2)));
        perpErr = [perpErr maxperpErr_ij];
      end
    end
end