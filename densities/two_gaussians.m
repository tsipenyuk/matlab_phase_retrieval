function [g_sol, sqrtI, X] = two_gaussians(nDim, varargin)
% two_gaussians - density corresponding to two uncorrelated gaussians
%
% Synopsis ([]s are optional)
%   [g_sol, sqrtI] = two_gaussians(nDim, [m1], [m2], [Sigma1], [Sigma2],...
%                                  [lb], [ub], [nPts])
%
% Description
%   Returns a density comprised of two gaussian densities with
%   means m1 and m2 and covariance matrices Sigma1 and Sigma2.
%   The densities are evaluated in the box specified by lower and 
%   upper bounds lb and ub; the box contains nPts points at which 
%   the density is approximated. 
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
%   (ndarray) [Sigma1]   Width of the first Gaussian
%             [ Sigma1 = [10] ] if nDim == 1
%             [ Sigma1 = [10 0; 0 10] ] if nDim == 2
%             [ Sigma1 = [100 0 0; 0 100 0; 0 0 100] ] if nDim == 3
%   (ndarray) [Sigma2]   Width of the second Gaussian
%             [ Sigma2 = [10] ] if nDim == 1
%             [ Sigma2 = [10 10] ] if nDim == 2
%             [ Sigma2 = [100 0 0; 0 100 0; 0 0 100] ] if nDim == 3
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
%      [g,~,x] = two_gaussians(2, 'm1', [0.6 0.3], 'm2', [-0.3 0], ...
%                          'Sigma1', [50 2; 2 10], ...
%                          'Sigma2', [25 -3; 18 10], ...
%                          'nPts', [200 200]);
%   contour(x{1},x{2},g);
%
%   %% 3D
%   [g,~,x] = two_gaussians(3, 'm1', [0.6 0.2 -0.2],...
%                           'm2', [-0.3 0 0.4], ...
%                          'Sigma1', [150 12 0; 12 110 0; 0 0 100], ...
%                          'Sigma2', [250 0 -13; 0 108 10; 13 203 200], ...
%                          'nPts', [70 70 70]);
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
        p.addParamValue('Sigma1', 10, @isscalar);
        p.addParamValue('Sigma2', 10, @isscalar);
        p.addParamValue('nPts', 100, @isscalar);
        p.addParamValue('lb', -1, @isscalar);
        p.addParamValue('ub', 1, @isscalar);
      case 2
        p.addParamValue('m1', [-0.5 -0.5], @isvector);
        p.addParamValue('m2', [0.5 0], @isvector);
        p.addParamValue('Sigma1', [10 0; 0 10], @ismatrix);
        p.addParamValue('Sigma2', [10 0; 0 10], @ismatrix);
        p.addParamValue('nPts', [100 100], @isvector);
        p.addParamValue('lb', [-1 -1], @isvector);
        p.addParamValue('ub', [1 1], @isvector);
      case 3
        p.addParamValue('m1', [-0.5 -0.5 0], @isvector);
        p.addParamValue('m2', [0.5 0 0.5], @isvector);
        p.addParamValue('Sigma1', [1e2 0 0; 0 1e2 0; 0 0 1e2], @ismatrix);
        p.addParamValue('Sigma2', [1e2 0 0; 0 1e2 0; 0 0 1e2], @ismatrix);
        p.addParamValue('nPts', [100 100 100], @isvector);
        p.addParamValue('lb', [-1 -1 -1], @isvector);
        p.addParamValue('ub', [1 1 1], @isvector);
    end
    p.parse(nDim, varargin{:});
    lb = p.Results.lb;
    ub = p.Results.ub;
    nPts = p.Results.nPts;
        
    x = arrayfun(@linspace, lb, ub, nPts, 'un', 0);
    % For higher dimensions, we remove singleton dimensions a lot
    s = @squeeze;
    % Switch-casing may be simplified with an extensive use of
    % arrayfun's, but the code becomes very unreadable. Since an
    % extension to nDim >= 4 is unlikely, this code is written as
    % it is.
    switch nDim
      case 1
        X1 = x{:};

        m1 = p.Results.m1;
        m2 = p.Results.m2;
        Sigma1 = p.Results.Sigma1;
        Sigma2 = p.Results.Sigma2;

        g_sol = exp(-(X1 - m1).^2 * Sigma1^2) + ...
                exp(-(X1 - m2).^2 * Sigma2^2);
        X = X1;
        
      case 2
        [X1, X2] = meshgrid(x{:});
        
        m1 = p.Results.m1;
        m2 = p.Results.m2;
        Sigma1 = p.Results.Sigma1;
        Sigma2 = p.Results.Sigma2;

        % Matrix multiplication -- calculate the exponent
        lng1 = zeros(size(X1));
        lng2 = zeros(size(X1));
        X = cat(3,X1,X2);
        
        for i=1:1:size(X1,1)
            for j=1:1:size(X1,2)
                lng1(i,j) = (s(X(i,j,:))' - m1) * ...
                            Sigma1 * (s(X(i,j,:)) - m1');
                lng2(i,j) = (s(X(i,j,:))' - m2) * ...
                            Sigma2 * (s(X(i,j,:)) - m2');
            end
        end
        g_sol = exp(-lng1) + exp(-lng2);
        X = {X1, X2};
        
      case 3
        [X1, X2, X3] = meshgrid(x{:});
        
        m1 = p.Results.m1;
        m2 = p.Results.m2;
        Sigma1 = p.Results.Sigma1;
        Sigma2 = p.Results.Sigma2;

        % Matrix multiplication -- calculate the exponent
        lng1 = zeros(size(X1));
        lng2 = zeros(size(X1));
        X = cat(4,X1,X2,X3);

        for i1=1:1:size(X1,1)
            disp(i1)
            for i2=1:1:size(X1,2)
                for i3=1:1:size(X1,3)
                    lng1(i1,i2,i3) = (s(X(i1,i2,i3,:))' - m1) * ...
                        Sigma1 * (s(X(i1,i2,i3,:)) - m1');
                    lng2(i1,i2,i3) = (s(X(i1,i2,i3,:))' - m2) * ...
                        Sigma2 * (s(X(i1,i2,i3,:)) - m2');
                end
            end
        end
        g_sol = exp(-lng1) + exp(-lng2);
        X = {X1, X2, X3};
    end
    
    sqrtI = abs(fftn(g_sol));    
end