function [gOut eOut] = compare_algorithms(g, sqrtI, nSteps, varargin)
% compare_algorithms - run multiple phase-ret. algoritms at once
%
% Synopsis ([]s are optional)
%   argout = compare_algorithms(g, sqrtI, nSteps, [varargin])
%
% Description
%   Performs nSteps of phase retrieval algorithm steps using
%   initial approximation g, data sqrtI (square root of intensity),
%   and various algorithms passed in varargin.
%
% Inputs ([]s are optional)
%   (ndarray) g      Current  approximation to the solution of the
%                    phase retrieval problem
%   (ndarray) sqrtI  Phase retrieval data (square root of the
%                    measured intensity)
%   (scalar)  nSteps Update parameter, default value (no explanation).
%   (...)     [varargin = {@er @hio_bauschke}]
%                    Cell array containing either:
%                      - handles to update functions, or
%                      - cells containing handles to update
%                      functions and additional arguments to update
%                      functions (see Examples). If an update
%                      function does not require any additional
%                      parameters, string 'None' must be provided
%                      as a parameter.
%
% Outputs
%   (cell array) gOut  Cell array approximated solutions of each algorithm;
%   (cell array) eOut  Cell array with energy descent information of each
%                         algorithm.
%
% Examples
%   %% 2D, two gaussians
%   [x1, x2] = ndgrid([-20:0.2:20], [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   A = abs(fftn(g_sol));
%   g_init = g;
%
%   % default comparisson
%   [gOut eOut] = compare_algorithms(g_init, A, 200);
%
%   % compare using two (or more) specific update algorithms:
%   [gOut eOut] = compare_algorithms(g_init, A, 200, @er, @hio_bauschke);
%
%   % Tell hio to use beta=0.6 and explicitely state pP as the
%   % projection in the object space:
%   [gOut eOut] = compare_algorithms(g_init, A, 200, {@er, 'None'}, ...
%                                            {@hio_bauschke, 0.6, @pP});
%   plot(gOut{1}) % Plot the approximation of er
%   plot(eOut{1}) % Plot the energy descent of er
%
% Requirements
%   If called with default parameters, the following functions are
%   required: 
%     pM (modulus projection)
%     pP (non-negative projection)
%     eM (modulus energy)
%     er (error-reduction update)
%     hio_bauschke (hybrid input-output update)
%  
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-08  First Edition
    
    % Parse input
    if nargin == 3
        disp('Comparing default algorithms: er and hio_bauschke...');
        algorithms = {@er, @hio_bauschke};
        params = {'None', 'None'};
    else
        if ~iscell(varargin{1}) % Parameters not specified
            algorithms = varargin;
            params = {};
            for iAlg = 1:1:length(algorithms)
                params{iAlg} = 'None';
            end
        else % Additional parameters specified
            algorithms = {};
            params = {};
            for iAlg = 1:1:length(varargin)
                algorithms{iAlg} = varargin{iAlg}{1};
                params{iAlg} = varargin{iAlg}{2};
            end
        end
    end

    % Initialize
    for iAlg = 1:1:length(algorithms)
        gOut{iAlg} = g;
        eNorm = eM(g, sqrtI);
        eOut{iAlg} = [1];
    end

    h = waitbar(0, '♚ ♛ ♜ ♝ ♞ ♟ ♔ ♕ ♖ ♗ ♘ ♙');
    % loop
    for iSteps = 1:1:nSteps
        h = waitbar(iSteps / nSteps);
        for iAlg = 1:1:length(algorithms)
            if params{iAlg} == 'None'
                [gOut{iAlg}, E] = ...
                    algorithms{iAlg}(gOut{iAlg}, sqrtI);
            else
                [gOut{iAlg}, E] = ...
                    algorithms{iAlg}(gOut{iAlg}, sqrtI, ...
                                     params{iAlg});
            end
            eOut{iAlg} = [eOut{iAlg} E / eNorm];
        end
    end
    
    argout = {gOut eOut};
    % Plot resulting energies 
    figure;
    ax = axes;
    m = ceil(sqrt(length(algorithms)));
    for iAlg = 1:1:length(algorithms)
        hold on;
        E = eOut{iAlg};
        plot(E);
        % Assuming the first algorithm is stabler, use it to fix
        % axes
        %axis([0 nSteps 0 1.1 * max(eOut{1}(:))])
    end
    set(ax ,'xscale','log');
    set(ax ,'yscale','log');
    % Convert algorithm names to plot legend string
    t = legend(cellstr(cellfun(@func2str, algorithms, 'un', 0)));
    % Allow underscores in function names
    set(t,'interpreter','none'); 
end