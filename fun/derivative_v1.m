function fout = derivative_v1(xin, Fin)
% NONREGULARDERIVATIVE one-dimensional periodic derivative
%
%   fout = nonregularDerivative(xin, Fin) is constructed as the inverse
%   function to the INDEFINITEINTEGRAL. It returns the derivative
%   corresponding to the discretized function at xin with
%   values Fin. The derivative is calculated using the formula
%   \begin{equation}
%   \label{eq1}
%   f_i + f_{i+1} = 2 (F_{i+1} - F_{i}) / (X_{i+1} - X_{i}).
%   \end{equation}
%   This relation does not set the first value of the derivative
%   $f_1$. To fix this value, we use the additional assumption that
%   the resulting derivative is periodic, i.e. that f_n = f_1.
%   This leads to the linear equation
%
%      / 1 1 0 0 ... 0 \   /      f_1        \
%      | 0 1 1 0 ... 0 |   |      f_2        |
%      | 0 0 1 1 ... 0 |   |      f_3        |
%      . 0 ...   ... 0 . x .      ...        . = RHS.
%      . 0 ...   ... 0 .   .      ...        .
%      | 0 ... 0 0 1 1 |   |      f_{n-2}    |
%      \ 1 0 0 ... 0 1 /   \      f_{n-1}    /
%    \______  __________/
%           \/
%           := A
%
%   and f_n = f_1. The matrix A is invertible if (n-1) is odd (n is even)
%   and may be computed directly to optimize performance. The
%   solution of \cref{eq1} remains valid if we add any random
%   number $a$ to $f_k$ and subtract $a$ from $f_{k+1}$ for all odd
%   $k$. We use this fact to regularize the resulting
%   derivative: we choose $a$ such that this derivative coincides
%   with the conventional two-point stencil derivative at xin(2).
%
%
%   xin = vector containing points at which the function is
%   discretized
%
%   Fin = vector containing function values at xin
%
%   fout = indefinite integral of fin discretized at xin
%
%   Arseniy Tsipenyuk, TUM M7
%   May 09, 2016
    if rem(length(Fin),2) == 1
        error('Length of each input vector must be even, not odd.')
    end
    
    if size(Fin,1) ~= 1
        Fin = Fin';
        xin = xin';
    end

    
    F1 = Fin(1:end-1);
    F2 = Fin(2:end);
    F1s = Fin(1:end-2);
    F3s = Fin(3:end);
    x1 = xin(1:end-1);
    x1s = xin(1:end-2);
    x3s = xin(3:end);
    x2 = xin(2:end);
    %x3 = circshift(xin, [0, 2]);
    
    RHS = 2 * (F2 - F1) ./ (x2 - x1);
    if size(RHS, 2) ~= 1
        RHS = RHS';
    end
    
    Mx = eye(length(RHS)) + circshift(eye(length(RHS)), [0,1]);
    size(Mx)
    
    f = inv(Mx) * RHS; % Vertical vector
    f = circshift(f', [0, -1]); % Horizontal vector
    
    % Regularize the result
    pm_ones = ones(size(f));
    pm_ones(2:2:end) = -1;
    %f(2)
    %a_new = sum((F3s - F1s) / (x3s - x1s) - f(2:end-1)) / (length(f)-2)
    f2 = (Fin(3) - Fin(1)) / (xin(3) - xin(1))
    a_new = f(2) - f2;
    f = f + a_new * pm_ones;
    fout = [f f(1)];
end
