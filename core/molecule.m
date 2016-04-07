classdef molecule
% MOLECULE describes a molecule (its density and its Fourier transform).
%
% Contains the following attributes:
%   obj.dimension - dimension of the problem (1, 2, or 3)
%   obj.nPoints   - array with size (1, obj.dimension), contains
%                   number of points in each dimension.
%   obj.xMin      - array with size (1, obj.dimension), contains
%                   lower bounds of the grid in the object space
%   obj.xMax      - --//--
%                   upper bounds of the grid in the object space
%   obj.kMin      - --//--
%                   lower bounds of the grid in the Fourier space
%   obj.kMax      - --//--
%                   upper bounds of the grid in the Fourier space
%   obj.xAxes     - cell array with size {1, obj.dimension},
%                   contains arrays of size (1,n), where n is a 
%                   corresponding number of points specified in 
%                   obj.nPoints. Elements of obj.xAxes are linear
%                   interpolations between obj.xMin and obj.xMax 
%                   in each dimension.
%   obj.kAxes     - cf. obj.xAxes; explanation applies mutatis
%                   mutandis.
%   obj.xMeshgrid - a full grid created using entries of xAxes;
%                   suited for evaluation of functions in the
%                   object space    
%   obj.kMeshgrid - --//--                    entries of kAxes
%                   --//--
%                   Fourier space    
%
% For examples, type 'doc mGrid' and see Constructor info.
% <NOT FINISHED>
    properties
        grid
        density
        ft
    end
    methods
        function obj=molecule(varargin)
        % MOLECULE Creates a molecule from density or fourier transform.
        %
        % Examples
        %
        % Creates a MOLECULE object (default: 1-dim, 200 pts., Gaussian):
        %    >> m = molecule()
        %    >> m.plot() % Look at the plot
        %
        % Creates a MOLECULE object (default: 1-dim, 200 pts., Gaussian):
        %    >> m = molecule()
        %    >> m.plot() % Look at the plot
        %
        % Evaluate function on a grid (2-dim. example)
        %    >> x = 0:2:16
        %    >> y = 0:3:50
        %    >> [X,Y] = meshgrid(x,y)
        %    >> f = sin(X + Y).^2
        %    >> m = molecule('density', f)
        %    >> surf(m.density)
        %
        % Additionally, define the grid
        %    >> x = 0:2:16
        %    >> y = 0:3:50
        %    >> [X,Y] = meshgrid(x,y)
        %    >> f = sin(X + Y).^2
        %    >> g = mGrid('xAxes',{x,y})
        %    >> m = molecule('density', f, 'grid', g)
        %    >> m.plot()
            
            p = inputParser;

            default_grid = mGrid(); % Default grid - 200 pts
                                    % between 0 and 2*pi
            default_density = exp( - (default_grid.xAxes{1}-pi) .^2);
            default_ft = fftn(default_density);
            
            addParameter(p,'grid',default_grid);
            addParameter(p,'density',default_density,@isnumeric);
            addParameter(p,'ft',default_ft,@isnumeric);

            parse(p,varargin{:});

            % Figure out which arguments were passed and which
            % construction method will be used
            str_varargin = varargin(1:2:size(varargin,2));
            if any(strmatch('density',str_varargin))
                obj.density = p.Results.density;
                obj.ft = fftn(obj.density);
                obj.grid = mGrid('nPoints', length(obj.density))
            elseif any(strmatch('ft',str_varargin))
                obj.ft = p.Results.ft;
                obj.density = real(ifftn(obj.ft));
                obj.grid = mGrid('nPoints', length(obj.ft))
            else
                disp(['Input density or FT were not specified. ' ...
                      'Proceeding with default density (Gaussian, ' ...
                      'width 1, mean pi, plotted between 0 and ' ...
                      '2*pi).']);
                obj.density = p.Results.density;
                obj.ft = p.Results.ft;
            end

            if any(strmatch('grid',str_varargin))
                obj.grid = p.Results.grid;
            else
                disp('(Warning: grid was not specified.)');
                if length(str_varargin) == 0
                    disp(['(Using the default grid, 200 pts between ' ...
                          '0 and 2*pi.)']);
                    obj.grid = p.Results.grid;
                else
                    disp(['(Using the default grid linearly spaced between ' ...
                          '0 and 2*pi.)']);
                end

            end
        end

        
        function set_density(obj, density)
            obj.density = density;
            obj.ft = fftn(obj.density);
        end

        
        function set_ft(obj, ft)
            obj.ft = ft;
            obj.density = real(ifftn(obj.ft));
        end

        
        function set_phase(obj, phase)
            obj.ft = abs(obj.ft) .* exp(1j * phase);
            obj.density = real(ifftn(obj.ft));
        end
        
        
        function set_random_phase(obj)
            random_phase = 2 * pi * rand(size(obj.ft));
            obj.set_phase(random_phase);
        end
        
        
        function set_zero_phase(obj, phase)
            obj.ft = abs(obj.ft);
            obj.density = real(ifftn(obj.ft));
        end


        function plot(obj, varargin)
        % PLOT plotting routine for molecule density and its FT
        %
        % Optional arguments: axes handles to axes where the
        % density and the FT modulus and phase will be plotted
        %
        % Examples
        % Simple plot
        %   >> m = molecule()
        %   >> m.plot()
        % 
        % Plot using specific axes
        %   >> figure();
        %   >> subplot(2,2,1)
        %   >> ax1 = gca
        %   >> subplot(2,2,4)
        %   >> ax2 = gca
        %   >> m.plot(ax1, ax2)
            if obj.grid.dimension == 1
                p = inputParser;
                
                addOptional(p,'ax1',0);
                addOptional(p,'ax2',0);
                
                parse(p,varargin{:});
                
                if length(varargin) == 2
                    ax1 = p.Results.ax1;
                    ax2 = p.Results.ax2;
                elseif length(varargin) == 0
                    fig = figure();
                    subplot(2,1,1)
                    ax1 = gca();
                    subplot(2,1,2)
                    ax2 = gca();
                else
                    disp(strcat(...
                        ['You must provide either two optional arguments ' ...
                         '(ax1, ax2, handles to axes where the ' ...
                         'density and FT  will be plotted), ' ...
                         'or none. ' ...
                         'Provided: '], ...
                        str(varargin), ', i.e. ', str(length(varargin)), ...
                        ' arguments.'))
                end
                
                hold on
                % 1st axis
                plot(ax1, obj.grid.xAxes{1}, obj.density)
                axis(ax1, [min(obj.grid.xAxes{1}) ...
                           max(obj.grid.xAxes{1}) ...
                           min(obj.density)...
                           max(obj.density)])
                ax1.XTick = [0 pi/2 pi  3*pi/2 2*pi];
                ax1.XTickLabel = {'0', '\pi/2', '\pi',  '3\pi/2', ...
                                  '2\pi'};
                title(ax1,'\fontsize{12}Density')

                
                
                % 2nd axis
                yyaxis left
                plot(ax2, obj.grid.kAxes{1}, abs(fftshift(obj.ft)))
                axis(ax2, [min(obj.grid.kAxes{1}) ...
                           max(obj.grid.kAxes{1}) ...
                           min(abs(fftshift(obj.ft))) ...
                           max(abs(fftshift(obj.ft)))])
                title(ax2,'\fontsize{12}FT: modulus and phase')
                
                % 3nd axis
                yyaxis right
                plot(ax2, obj.grid.kAxes{1}, ...
                     angle(fftshift(obj.ft)))
                axis(ax2, [min(obj.grid.kAxes{1}) ...
                           max(obj.grid.kAxes{1}) ...
                           -pi ...
                           pi])

                hold off
                
            elseif obj.grid.dimension == 2
                surf(obj.grid.xMeshgrid{1}, obj.grid.xMeshgrid{2}, ...
                     obj.density)
            else
                disp(['Currently, only 1-dimensional plotting method ' ...
                      'is supported'])
            end
        end
        
    end
end