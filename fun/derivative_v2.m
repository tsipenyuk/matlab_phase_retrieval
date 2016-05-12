function fout = derivative_v2(xin, Fin)
% DERIVATIVE one-dimensional periodic derivative 
%
%   fout = derivative(xin, Fin) is constructed as the inverse
%   function to the INDEFINITEINTEGRAL. It returns the derivative
%   corresponding to the discretized function at xin with
%   values Fin. The derivative is calculated using the formula
%   \begin{equation}
%   \label{eq1}
%   f_i + f_{i+1} = \underbrace{2 (F_{i+1} - F_{i}) / (X_{i+1} -
%                    X_{i})}_{RHS_i}.
%   \end{equation}
%   This relation does not set the first value of the derivative
%   $f_1$. To fix this value, we use the additional assumption that
%   the resulting derivative is periodic, i.e. that f_n = f_1. To
%   ensure the regularity of the result, we use the corresponding
%   second-order equation, namely
%   \begin{equation}
%   \label{eq1}
%   f_i + 2*f_{i+1} + f_{i + 2} = RHS_i + RHS_{i+1}.
%   \end{equation}
%   This leads to the linear equation
%
%      / 1 2 1 0 ... 0 \   /      f_1        \
%      | 0 1 2 1 ... 0 |   |      f_2        |
%      | 0 0 1 2 ... 0 |   |      f_3        |
%      . 0 ...   ... 0 . x .      ...        . = RHS_i + RHS_{i+1}.
%      . 0 ...   1 2 1 .   .      ...        .
%      | 1 ... 0 0 1 2 |   |      f_{n-2}    |
%      \ 2 1 0 ... 0 1 /   \      f_{n-1}    /
%    \______  __________/
%           \/
%           := A
%
%   and f_n = f_1. The matrix A is invertible if (n-1) is odd (n is even)
%   and may be computed directly to optimize performance.
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
%   May 11, 2016
    if rem(length(Fin),2) ~= 0
        error('Length of each input vector must be even, not odd.')
    end

    F1 = Fin(1:end-1);
    F2 = Fin(2:end);
    x1 = xin(1:end-1);
    x2 = xin(2:end);

    RHS = 2 * (F2 - F1) ./ (x2 - x1);
    if size(RHS, 2) ~= 1
        RHS = RHS';
    end
    RHS = RHS + circshift(RHS, [1,0]);

    Mx = eye(length(RHS)) + 2 * circshift(eye(length(RHS)), [0,1]) + 1 * circshift(eye(length(RHS)), [0,2]);
    
    f = inv(Mx) * RHS; % Vertical vector
    fout = [f' f(1,1)]; % Horizontal vector
    
    % Regularize the result
    %a = (fout(2) - fout(1)) / 2;
    %pm_ones = ones(size(fout));
    %pm_ones(2:2:end) = -1;
    %fout = fout + a * pm_ones; 
end