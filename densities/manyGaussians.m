function [gSol, sqrtI, x, s] = manyGaussians(nDim, nGauss, varargin)
% manyGaussians - density corresponding to nGauss uncorrelated gaussians
%
% Synopsis ([]s are optional)
%   [gSol, sqrtI, X, s] = many_gaussians(nDim, nGauss, 
%                                        [lb], [ub], [nPts])
%
% Description
%   Returns a density 
%
% Inputs ([]s are optional)
%   (integer) nDim         Number of dimensions, set to 1, 2, or 3
%   (integer) nGauss       Number of gaussians
%   (integer) [nPts = 200] Default number of points in the
%                          discretization grid
%   (integer) [randomSeed] State of the random seed generator,
%                          left unchanged by default
%   (float)   [squeezeBox = 0.5]   
%                          To ensure that the density at the
%                          boundaries is zero, means of gaussians
%                          are generated inside a subset (box) of
%                          the total discretized domain. The
%                          distance between the box and the domain
%                          boundary is defined by squeezeBox.
%   (float)   [varMean]
%                          Mean value of the gaussian variance. By
%                          default chosen to be 1/3 * (expected average
%                          distance between gaussians).
%   (float)   [varVar]
%                          Variance value of the gaussian variance. By
%                          default chosen to be 1/3 * (expected average
%                          distance between gaussians).
% Outputs
%   (ndarray) gSol         Solution density
%   (ndarray) sqrtI        Its Fourier transform modulus
%   (ndarray) x            Cell array of meshgrid coordinates at
%                          which gSol is evaluated 
%   (integer) s            Random seed from which the density was generated
%
% Examples
%   %% 1D
%   [g,~,x,s] = manyGaussians(1, 5);
%   plot(x{1},g);
%
%   %% 2D
%   [g,~,x,s] = manyGaussians(2,5);
%   contour(x{1},x{2},g);
%   %% Repeat the same result
%   [g,~,x,s] = manyGaussians(2, 5, 'randomSeed', s);
%   contour(x{1},x{2},g);
%   %% Control gaussian variance parameters   
%   [g,~,x,s] = manyGaussians(2,5, 'varMean', 0.2, 'varVar', 0.1);
%   contour(x{1},x{2},g);
%
%
%   %% 3D
%   [g,~,x,s] = manyGaussians(2,5);
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
    squeezeVarFactor = 8;
    
    % Mean of the variance of the gaussians
    varMeanDefault =  avgDistance / squeezeVarFactor;
    % Variance of the variance of the gaussians
    varVarDefault = varMeanDefault / 2;

    
    p = inputParser;
    p.addRequired('nDim', @(x) any(x==[1 2 3]));
    p.addRequired('nGauss', @isfloat);
    p.addParamValue('nPts', nPtsDefault, @isvector);
    p.addParamValue('randomSeed', 0, @isfloat);
    p.addParamValue('squeezeBox', 0.5, @isfloat);
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
    m = (1 - pR.squeezeBox) * (2 * rand([nGauss nDim]) - 1)
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