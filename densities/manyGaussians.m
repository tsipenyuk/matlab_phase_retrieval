function [gSol, sqrtI, X, s] = manyGaussians(nDim, nGauss, varargin)
% manyGaussians - density corresponding to nGauss uncorrelated gaussians
%
% Synopsis ([]s are optional)
%   [gSol, sqrtI] = two_gaussians(nDim, nGauss, 
%                                  [lb], [ub], [nPts])
%
% Description
%   Returns a density 
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
%   (ndarray)     gSol  Molecule density
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
%   2016-06-28  First Edition

    % Default number of points in the discretization grid
    nPtsDefault = 200 * ones([1, nDim]);
    
    % Average distance between nGauss points placed in a [-1 1] box
    avgDistance = 2 * sqrt(nDim) / nGauss;
    squeezeVarFactor = 3;
    
    % Mean of the variance of the gaussians
    varMeanDefault =  avgDistance / squeezeVarFactor;
    % Variance of the variance of the gaussians
    varVarDefault = avgDistance / squeezeVarFactor;

    
    p = inputParser;
    p.addRequired('nDim', @(x) any(x==[1 2 3]));
    p.addRequired('nGauss', @isfloat);
    p.addParamValue('nPts', nPtsDefault, @isvector);
    p.addParamValue('randomSeed', 0, @isfloat);
    p.addParamValue('squeezeBox', 3, @isfloat);
    p.addParamValue('varMean', varMeanDefault, @isfloat);
    p.addParamValue('varVar', varVarDefault, @isfloat);
    p.parse(nDim, nGauss, varargin{:});
    pR = p.Results;
    
    % Set random seed, if not specified
    if ismember('randomSeed', p.UsingDefaults)
        rng('shuffle');
        s = randi(450638601); % Any large number to get enough states
        rng(s, 'twister');
        disp(['Random seed set to ' num2str(s) '.']);
    else
        s = pR.randomSeed;
        rng(s, 'twister');
    end
    
    % Define upper/lower bounds of the box where the density is evaluated
    switch nDim
      case 1
        lb = -1;
        ub = 1;
      case 2
        lb = [-1 -1];
        ub = [ 1  1];
      case 3
        lb = [-1 -1 -1];
        ub = [ 1  1  1];
    end
    
    % Density is evaluated at X
    x = arrayfun(@linspace, lb, ub, pR.nPts, 'un', 0);
    % Generate random means/variances of gaussians
    m = (1 - pR.squeezeBox * pR.varMean) * (2 * rand([nGauss nDim]) - 1)
    v = abs(normrnd(pR.varMean, pR.varVar, [nGauss nDim]))
    

    % Function that evaluates the density
    function gSol = density(x, m, v)
        switch pR.nDim
          case 1
            X = x{1};
            gSol = zeros(size(X));
            for iGauss = 1:pR.nGauss
                gSol = gSol + ...
                       exp( - (X - m(iGauss)).^2 / v(iGauss).^2);
            end
          case 2
            [X Y] = ndgrid(x{1}, x{2});
            gSol = zeros(size(X));
            for iGauss = 1:pR.nGauss
                gSol = gSol + ...
                       exp( - (X - m(iGauss,1)).^2 / v(iGauss,1).^2 + ...
                            - (Y - m(iGauss,2)).^2 / v(iGauss,2).^2);
            end
          case 3
            [X, Y, Z] = ndgrid(x{1}, x{2}, x{3});
            gSol = zeros(size(X));
            for iGauss = 1:pR.nGauss
                gSol = gSol + ...
                       exp( - (X - m(iGauss,1)).^2 / v(iGauss,1).^2 + ...
                            - (Y - m(iGauss,2)).^2 / v(iGauss,2).^2 + ...
                            - (Z - m(iGauss,3)).^2 / v(iGauss,3).^2);
            end
        end
    end
    
    gSol = density(x, m, v);
    sqrtI = abs(fftn(gSol));
end