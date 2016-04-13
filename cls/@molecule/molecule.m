classdef molecule
% MOLECULE describes a molecule (its density and its Fourier transform).
%
% Contains the following attributes:
%   obj.grid      - class mGrid, stores dimension and grids
%   obj.density   - array storing the dencity of the molecule
%   obj.ft      - array storing the Fourier transform of the molecule
%
% For examples, type 'doc molecule' and see Constructor info.
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
        % When the instance of the class is created, one must provide
        % either the density or the FT; otherwise, the default
        % 1-dim. density is used. If the density is
        % one-dimensional, a 'grid' object is created
        % automatically.
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
                obj.grid = mGrid('nPoints', length(obj.density));
            elseif any(strmatch('ft',str_varargin))
                obj.ft = p.Results.ft;
                obj.density = real(ifftn(obj.ft));
                obj.grid = mGrid('nPoints', length(obj.ft));
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
                if isvector(obj.density)
                    if length(str_varargin) == 0
                        disp(['(Using the default grid, 200 pts between ' ...
                              '0 and 2*pi.)']);
                        obj.grid = p.Results.grid;
                    else
                        disp(['(Using the default grid linearly spaced '...
                              'between 0 and 2*pi.)']);
                        obj.grid = mGrid('nPoints', length(obj.density))
                    end
                else 
                    disp(['(Creating a default grid; in object space, '...
                           'box with linearly spaced stencils between '...
                           '0 and 2pi.)'])
                    dim = length(size(obj.density));
                    % For two axes x1 and x2 with lengths i and j,
                    % meshgrid yields
                    % result with size [j  i]. fliplr() accounts
                    % for that dimension mixture
                    custom_size = fliplr(size(obj.density));
                    custom_nPoints = custom_size;
                    custom_xMin = zeros(1, dim);
                    custom_xMax = 2 * pi * ones(1, dim);
                    custom_kMin = ones(1, dim);
                    custom_kMax = custom_size;
                    grid = mGrid('nPoints', custom_nPoints, ...
                                 'xMin', custom_xMin, 'xMax', custom_xMax, ...
                                 'kMin', custom_kMin, 'kMax', ...
                                 custom_kMax); 
                    obj.grid = grid;
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


        % Longer functions are specified in their own files
        plot(obj, varargin);
        contour(obj, varargin);
    end
end