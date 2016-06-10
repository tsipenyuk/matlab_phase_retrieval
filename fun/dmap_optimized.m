function [g_new, error] = dmap_optimized(g, sqrtI, beta, pObj, varargin)
% hio - Optimized difference map algorithm [Elser]
%
% Synopsis ([]s are optional)
%   [g_new, error] = dmap_optimized(g, sqrtI, [beta], [pObj], [varargin])
%
% Description
%   Performs a phase retrieval update step and calculates the error
%   (energy) corresponding to the optimized difference map update,
%   cf. [Elser], eq. 18.
%
% Inputs ([]s are optional)
%   (ndarray) g      Current  approximation to the solution of the
%                    phase retrieval problem
%   (ndarray) sqrtI      Phase retrieval data (square root of the
%                    measured intensity)
%   (scalar)  [beta = 1]
%                    Update parameter, default value reproduces the
%                    hybrid input-output update (Bauschke' variant)
%   (func)    [pObj = @pP]
%                    handle to the projection onto the object space
%                    (non-negativity, support, etc.). pObj must
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
%   sqrtI = abs(fftn(g_sol));
%   g_new = zeros(size(g));
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap_optimized(g_new, sqrtI);
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
%   sqrtI = abs(fftn(g_sol));
%   g_new = zeros(size(g));
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap_optimized(g_new, sqrtI);
%       E = [E error];
%   end
%
% See also
%   er
%   bio
%   hio_fienup
%   dmap
%   
% Requirements
%   pM (modulus projection)
%   pP (non-negative projection)
%
% References
%   V. Elser, “Phase retrieval by iterated projections,” 
%   J. Opt. Soc. Am. A., vol. 20, pp. 40–55, 2003.
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-10  First Edition
    if nargin == 2
        beta = 1;
        pObj = @pP;
    end
    
    if nargin == 3
        pObj = @pP;
    end
    
    % Calculate the update
    if nargin <= 4
        F2 = (1 + 1. / beta) * pM(g, sqrtI) - 1. / beta * g;
        F1 = (1 - 1. / beta) * pObj(g) + 1. / beta * g;
        D = pObj(F2) - pM(F1, sqrtI);
        g_new = g + beta * D;
        error = eM(pObj(g_new), sqrtI);
    else
        F2 = (1  + 1. / beta) * pM(g, sqrtI) - 1. / beta * g;
        F1 = (1 - 1. / beta) * pObj(g, varargin) + 1. / beta * g;
        D = pObj(F2, varargin) - pM(F1, sqrtI);
        g_new = g + beta * D;
        error = eM(pObj(g_new, varargin), sqrtI);
    end
    g_new = pM(g_new, sqrtI);
end