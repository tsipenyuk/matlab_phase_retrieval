function [g_new, error] = bio(g, A, beta)
% bio - Basic Input-Output algorithm
%
% Synopsis ([]s are optional)
%   [g_new, error] = bio(g, A, [beta], [pObj], [varargin])
%
% Description
%   Performs a phase retrieval update step and calculates the error
%   (energy) corresponding to the update.
%
% Inputs ([]s are optional)
%   (ndarray) g      Current  approximation to the solution of the
%                    phase retrieval problem
%   (ndarray) A      Phase retrieval data (square root of the
%                    measured intensity)
%   (scalar)  [beta = 1.25]
%                    Update parameter, default value from [Fienup]
%   (func)    [pObj = @pP]
%                    handle to the projection onto the object space
%                    (non-negativity, atomicity, etc.). pObj must
%                    take g as the first argument, may contain
%                    optional arguments specified in varargin.
%   (...)     [varargin]
%                    optional arguments passed to pObj --- if
%                    submitted, pObj is  called as pObj(g, varargin)
%
% Outputs ([]s are optional)
%   (ndarray) g_new  Updated approximation to the solution of the 
%                    phase retrieval problem
%   (scalar)  error  Error (energy) corresponding to the updated 
%                    approximation
%
% Examples
%   %% 1D, two gaussians
%   x = [-20:0.2:20];
%   g_sol = exp(-x.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   A = abs(fftn(g_sol));
%   g_new = A;
%   E = [];
%   for i=1:1:200
%       % Stabilizing ER step, cf. [Fienup], p. 2765
%       [g_new, error] = er(g_new, A); 
%       E = [E error];
%       [g_new, error] = bio(g_new, A);
%       E = [E error];
%   end
%   plot(E);
%   plot(g_new);
%
%   %% 2D, two gaussians
%   [x1, x2] = ndgrid([-20:0.2:20], [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   A = abs(fftn(g_sol));
%   g_new = A;
%   E = [];
%   for i=1:1:200
%       % Stabilizing ER step, cf. [Fienup], p. 2765
%       [g_new, error] = er(g_new, A); 
%       E = [E error];
%       [g_new, error] = bio(g_new, A);
%       E = [E error];
%   end
%
% References
%   J. R. Fienup, “Phase retrieval algorithms: a comparison,” 
%       Applied Optics, vol. 21, pp. 2758–2769, 1982.
%   doc/phase_retrieval_algorithms.pdf
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase Retrieval Sandbox root folder.
%
% Changes
%   2016-06-01  First Edition
%   2016-06-07  Added pObj support
    if nargin == 2
        beta = 1.25;
        pObj = @pP;
    end
    
    if nargin == 3
        pObj = @pP;
    end
    
    % Calculate the update
    pM_g = pM(g, A);
    if nargin <= 4
        g_new = g - beta * pM_g + beta * pObj(pM_g);
    else
        g_new = g - beta * pM_g + beta * pObj(pM_g, varargin);
    end
    error = eM(g_new, A);
end