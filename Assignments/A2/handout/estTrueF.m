function [F0] = estTrueF(dLeft, dRight, Rleft, Rright, MintLeft, MintRight)
    % Estimate True Fundamental Matrix, F0
    % Tvec: dLeft - dRight
    % T: turn cross product into matrix multiplication
    % E: Essential Matrix
    % F0: True Fundamental Matrix (inv(MintLeft)' * E * inv(MintRight))
    
    Tvec = dLeft - dRight;
    T = zeros(3,3);
    T(1) = 0;
    T(2) = Tvec(3);
    T(3) = Tvec(2);
    T(4) = Tvec(3);
    T(5) = 0;
    T(6) = -Tvec(1);
    T(7) = -Tvec(2);
    T(8) = Tvec(1);
    T(9) = 0;
    E = Rleft* T * Rright';                       % Essential matrix
    F0 = inv(MintLeft)' * E * inv(MintRight);     % Ground-truth F matrix
end