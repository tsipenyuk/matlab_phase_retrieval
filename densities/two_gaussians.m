function [g_sol, sqrtI, X] = two_gaussians(nDim, varargin)
% two_gaussians - density corresponding to two uncorrelated gaussians
%
% Synopsis ([]s are optional)
%   [g_sol, sqrtI] = two_gaussians(nDim, [m1], [m2], [w1], [w2],...
%                                  [lb], [ub], [nPts])
%
% Description
%   Returns a density comprised of two gaussian densities with
%   means m1 and m2 and widths w1 and w2, assuming that there is no
%   correlation between the variables of the gaussians. The densities
%   are evaluated in the box specified by lower and upper bounds lb
%   and ub; the box contains nPts points at which the density is
%   approximated. 
%
% Inputs ([]s are optional)
%   (scalar)  nDim   Number of dimensions, set to 1, 2, or 3
%   (ndarray) [m1]   Mean of the first Gaussian
%             [ m1 = [-0.5] ] if nDim == 1
%             [ m1 = [-0.5 -0.5] ] if nDim == 2
%             [ m1 = [-0.5 -0.5 0] ] if nDim == 3
%   (ndarray) [m2]   Mean of the second Gaussian
%             [ m2 = [0.5] ] if nDim == 1
%             [ m2 = [0.5 0] ] if nDim == 2
%             [ m2 = [0.5 0 0.5] ] if nDim == 3
%   (ndarray) [w1]   Width of the first Gaussian
%             [ w1 = [0.1] ] if nDim == 1
%             [ w1 = [0.1 0.1] ] if nDim == 2
%             [ w1 = [0.1 0.1 0.1] ] if nDim == 3
%   (ndarray) [w2]   Width of the second Gaussian
%             [ w2 = [0.1] ] if nDim == 1
%             [ w2 = [0.1 0.1] ] if nDim == 2
%             [ w2 = [0.1 0.1 0.1] ] if nDim == 3
%   (ndarray) [nPts] Number of points in each dimension
%             [ nPts = [100] ] if nDim == 1
%             [ nPts = [100 100] ] if nDim == 2
%             [ nPts = [100 100 100] ] if nDim == 3
%   (ndarray) [lb]   Lower bound of the coordinates
%             [ lb = [-1] ] if nDim == 1
%             [ lb = [-1 -1] ] if nDim == 2
%             [ lb = [-1 -1 -1] ] if nDim == 3
%   (ndarray) [ub]   Upper bound of the coordinates
%             [ ub = [1] ] if nDim == 1
%             [ ub = [1 1] ] if nDim == 2
%             [ ub = [1 1 1] ] if nDim == 3
%
% Outputs ([]s are optional)
%   (ndarray)     g_sol  Molecule density
%   (ndarray)     sqrtI  Fourier modulus of the molecule density
%   (cell array)  X      Array of gridvectors xgv (ygv, zgv) at which
%                        the density is evaluated. May be used to
%                        calculate meshgrid(X{1}, X{2}) or similar.
%
% Examples
%   %% 1D
%   [g,~,x] = two_gaussians(1);
%   plot(x,g);
%
%   %% 2D
%   [g,~,x] = two_gaussians(2);
%   contour(x{1},x{2},g);
%
%   %% 2D with specified means, widths, and number of points
%  [g,~,x] = two_gaussians(2, 'm1', [0.8 0], 'm2', [-0.3 0], ...
%                         'w1', [0.05 0.05], 'w2', [0.1 0.25], ...
%                         'nPts', [500 500]);
%
%   %% 3D
%   [g,~,x] = two_gaussians(3);
%   plot_isosurface(x{1},x{2},x{3},g); % requires ../fun
%
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-09  First Edition

    p = inputParser;
    p.addRequired('nDim', @(x) any(x==[1 2 3]));
    
    % Fill in unset optional values.
    switch nDim
      case 1
        p.addParamValue('m1', -0.5, @isscalar);
        p.addParamValue('m2', 0.5, @isscalar);
        p.addParamValue('w1', 0.1, @isscalar);
        p.addParamValue('w2', 0.1, @isscalar);
        p.addParamValue('nPts', 100, @isscalar);
        p.addParamValue('lb', -1, @isscalar);
        p.addParamValue('ub', 1, @isscalar);
      case 2
        p.addParamValue('m1', [-0.5 -0.5], @isvector);
        p.addParamValue('m2', [0.5 0], @isvector);
        p.addParamValue('w1', [0.1 0.1], @isvector);
        p.addParamValue('w2', [0.1 0.1], @isvector);
        p.addParamValue('nPts', [100 100], @isvector);
        p.addParamValue('lb', [-1 -1], @isvector);
        p.addParamValue('ub', [1 1], @isvector);
      case 3
        p.addParamValue('m1', [-0.5 -0.5 0], @isvector);
        p.addParamValue('m2', [0.5 0 0.5], @isvector);
        p.addParamValue('w1', [0.1 0.1 0.1], @isvector);
        p.addParamValue('w2', [0.1 0.1 0.1], @isvector);
        p.addParamValue('nPts', [100 100 100], @isvector);
        p.addParamValue('lb', [-1 -1 -1], @isvector);
        p.addParamValue('ub', [1 1 1], @isvector);
    end
    p.parse(nDim, varargin{:});
    lb = p.Results.lb;
    ub = p.Results.ub;
    nPts = p.Results.nPts;
        
    x = arrayfun(@linspace, lb, ub, nPts, 'un', 0);
    % Switch-casing may be simplified with an extensive use of
    % arrayfun's, but the code becomes very unreadable. Since an
    % extension to nDim >= 4 is unlikely, this code is written as
    % it is.
    switch nDim
      case 1
        X1 = x{:};

        m1 = p.Results.m1;
        m2 = p.Results.m2;
        w1 = p.Results.w1;
        w2 = p.Results.w2;

        g_sol = exp(-(X1 - m1).^2 / w1^2) + ...
                exp(-(X1 - m2).^2 / w2^2);
        X = X1;
        
      case 2
        [X1, X2] = meshgrid(x{:});
        m1 = p.Results.m1;
        m2 = p.Results.m2;
        w1 = p.Results.w1;
        w2 = p.Results.w2;

        g_sol = exp(-(X1 - m1(1)).^2 / w1(1)^2 ...
                    -(X2 - m1(2)).^2 / w1(2)^2) + ...
                exp(-(X1 - m2(1)).^2 / w2(1)^2 ...
                    -(X2 - m2(2)).^2 / w2(2)^2);
        X = {X1, X2};
        
      case 3
        [X1, X2, X3] = meshgrid(x{:});
        
        m1 = p.Results.m1;
        m2 = p.Results.m2;
        w1 = p.Results.w1;
        w2 = p.Results.w2;
        
        g_sol = exp(-(X1 - m1(1)).^2 / w1(1)^2 ...
                    -(X2 - m1(2)).^2 / w1(2)^2 ...
                    -(X3 - m1(3)).^2 / w1(3)^2) + ...
                exp(-(X1 - m2(1)).^2 / w2(1)^2 ...
                    -(X2 - m2(2)).^2 / w2(2)^2 ...
                    -(X3 - m2(3)).^2 / w2(3)^2);
        X = {X1, X2, X3};
    end
    
    sqrtI = abs(fftn(g_sol));    
end