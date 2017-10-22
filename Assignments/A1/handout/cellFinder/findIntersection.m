function [intersection] = findIntersection(x1, n1, x2, n2)
    % Find intersection between a line x1 + p*n1 and
    % x2 + q*n2
    % x1: 2x1 array, denotes a 2D point
    % n1: 2x1 array, denotes a 2D normal vector
    % x2: 2x1 array, denotes a 2D point
    % n2: 2x1 array, denotes a 2D normal vector
    % line1 = x1 + p*n1
    % line2 = x2 + q*n2
    % Find p such that line1 = line2
    
    p = ((x2(1) - x1(1)) * n2(2) - (x2(2) - x1(2)) * n2(1)) ...
                        / (n1(1) * n2(2) - n1(2) * n2(1));
    if abs(n1(1) * n2(2) - n1(2) * n2(1)) < 1e-2
        isparalell = true;
    else
        isparalell = false;
    end
    
    if isparalell                    % If the two lines are paralell: det=0
        intersection = (x1 + x2) / 2;% return middle of the two points
    else
        intersection = x1 + p * n1;  % return intersection point
    end
end