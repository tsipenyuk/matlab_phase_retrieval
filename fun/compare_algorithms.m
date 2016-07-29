function [gOut eOut M] = compare_algorithms(g, sqrtI, nSteps, varargin)
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
%   %% 2D, five gaussians
%   [g,sqrtI,x,~] = ...
%   manyGaussians(2,15,'randomSeed',283362422,'varMean', 0.05,...
%   'varVar', 0.001); contour(x{1},x{2},g);
%   sqrtIpadded = ones(size(sqrtI));
%   sqrtIpadded(~isnan(sqrtI)) = sqrtI(~isnan(sqrtI));
%   gInit = pP(real(ifftn(sqrtIpadded.*sqrtIpadded)));
%    [gOut eOut] = compare_algorithms(gInit, sqrtI, 1000, @er, @hio_bauschke, @hio_simplified, @hio_fienup);
%
%   % default comparisson
%   [gOut eOut] = compare_algorithms(gInit, A, 200);
%
%   % compare using two (or more) specific update algorithms:
%   [gOut eOut] = compare_algorithms(gInit, A, 200, @er, @hio_bauschke);
%
%   % Tell hio to use beta=0.6 and explicitely state pP as the
%   % projection in the object space:
%   [gOut eOut] = compare_algorithms(gInit, A, 200, @er, ...
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
    
    DEBUG = 1;
    
    % Parse input
    if nargin == 3
        disp('Comparing default algorithms: er and hio_bauschke...');
        algorithms = {@er, @hio_bauschke};
        params = {'None', 'None'};
        areParamSpecified = {false false};
    else
        algorithms = {};
        params = {};
        areParamSpecified = {};
        for iVar = 1:1:length(varargin)
            if ~iscell(varargin{iVar}) % Parameters not specified
                areParamSpecified{iVar} = false;
                algorithms{iVar} = varargin{iVar};
                params{iVar} = 'None';
            else % Additional parameters specified
                areParamSpecified{iVar} = true;
                algorithms{iVar} = varargin{iVar}{1};
                params{iVar} = {varargin{iVar}{2:end}};
        end
    end

    % Initialize
    for iAlg = 1:1:length(algorithms)
        gOut{iAlg} = g;
        eNorm = eM(g, sqrtI);
        eOut{iAlg} = [1];
    end

    h = waitbar(0, '♚ ♛ ♜ ♝ ♞ ♟ ♔ ♕ ♖ ♗ ♘ ♙');
    h2 = figure;
    ax2 = axes;
    
    %DEBUG
    if DEBUG == 1
        sparseFactor = 10;
        nFrames = fix(nSteps/sparseFactor);
        
        global gSol; %global a; global b;
                     %F(nSteps) = struct('cdata',[],'colormap',[]);
        hDebug = figure;
        axDebug = gca();
        set(hDebug, 'Visible', 'off');
        
        M = moviein(nFrames);
        set(gca,'NextPlot','replacechildren');
    else 
        M = {}
    end
    % loop
    for iSteps = 1:1:nSteps
        h = waitbar(iSteps / nSteps);
        for iAlg = 1:1:length(algorithms)
            if areParamSpecified{iAlg} == false
                [gOut{iAlg}, E] = ...
                    algorithms{iAlg}(gOut{iAlg}, sqrtI);
            else
                [gOut{iAlg}, E] = ...
                    algorithms{iAlg}(gOut{iAlg}, sqrtI, ...
                                     params{iAlg}{:});
            end
            eOut{iAlg} = [eOut{iAlg} E / eNorm];
        end
        
        if rem(iSteps, 200) == 0
            cla(ax2);
            m = ceil(sqrt(length(algorithms)));
            for iAlg = 1:1:length(algorithms)
                hold on;
                E = eOut{iAlg};
                plot(ax2, E);
            end
            set(ax2 ,'xscale','log');
            set(ax2 ,'yscale','log');
            % Convert algorithm names to plot legend string
            t = legend(cellstr(cellfun(@func2str, algorithms, 'un', 0)));
            % Allow underscores in function namesn
            set(t,'interpreter','none'); 
        end
        
        % DEBUG
        if DEBUG == 1
            if rem(iSteps,10) == 0
                iFrame = fix(iSteps/sparseFactor)
                % Requires global variable 'g' describing the solution
                % of the problem
                figure(hDebug);
                set(hDebug,'visible','off');
                subplot(3,1,1); 
                hold on;
                plot(gSol);
                plot(gOut{1});
                legend('solution', 'er approximation')
                hold off;
                subplot(3,1,2);
                hold on;
                plot(unwrap(fftshift(angle(fftn(gSol)))));
                yyaxis right;
                plot(unwrap(fftshift(angle(fftn(gOut{1}))))); 
                hold off; 
                legend('solution phase', 'er phase');
                subplot(3,1,3);
                plot(fftshift(abs(fftn(gSol))));

                saveas(hDebug, strcat('output/', num2str(iSteps), ...
                                  '.png'))
                M(:,iFrame) = getframe;
                clf(hDebug,'reset');
            end
        end
        % DEBUG
    end
    
    if DEBUG == 1
        close(hDebug)
        %movie2avi(M, 'output/movie.avi')
    end

    close(h);
    
    argout = {gOut eOut};
    % Plot resulting energies 
    %figure;
    %ax = axes;
    %
    %m = ceil(sqrt(length(algorithms)));
    %for iAlg = 1:1:length(algorithms)
    %    hold on;
    %    E = eOut{iAlg};
    %    plot(E);
    %end
    %set(ax ,'xscale','log');
    %set(ax ,'yscale','log');
    %% Convert algorithm names to plot legend string
    %t = legend(cellstr(cellfun(@func2str, algorithms, 'un', 0)));
    %% Allow underscores in function namesn
    %set(t,'interpreter','none'); 
end