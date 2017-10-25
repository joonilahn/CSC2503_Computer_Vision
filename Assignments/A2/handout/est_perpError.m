function [max_perpErr] = est_perpError(x, y, F0, F1, nPoints)
% Calculate the maximum perpendicular error
    x = linspace(-150, 150, nPoints);
    y = linspace(-100, 100, nPoints);
    [X, Y] = meshgrid(x, y);
    axis([-150 150 -100 100]); axis xy;
    ax = axis;
    cropBox = [ax(1) ax(3) ax(2) ax(4)];

    perpErr = [];
    for i = 1:nPoints
        for j = 1:nPoints
          lF0 = [x(i); y(j); 1]' * F0';
          lF1 = [x(i); y(j); 1]' * F1';
          ep0 = cropLineInBox(lF0(1:2), lF0(3), cropBox);
          if sum(sum(isnan(ep0))) > 0
              break;
          else
            maxperpErr_ij = max(abs([ep0 ones(2,1)]*lF1') / norm(lF1(1:2)));
            perpErr = [perpErr maxperpErr_ij];
          end
        end
    end
    max_perpErr = max(perpErr);
end