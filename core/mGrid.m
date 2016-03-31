classdef mGrid
% MGRID describes the coordinates on which the density of a molecule
% and its transform are evaluated.
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
%
% For examples, type 'doc mGrid' and see Constructor info.
    properties
        dimension
        nPoints
        xMin
        xMax
        kMin
        kMax
        xAxes
        kAxes
    end
    methods
        function obj=mGrid(varargin)
        % MGRID Creates a grid on which density and its FT are evaluated.
        %
        % There are two distinct construction methods implicitely
        % specified in this constructor. In the first method, one
        % specifies the number of points, the lowest, and the
        % highest values of the grid. In the second method, one
        % directly inputs the grids.
        %
        % Examples
        % First method examples:
        %
        % Creates a MGRID object (default: 1-dim, 200 pts.):
        %    >> g = mGrid()
        %    >> g.xAxes{1,:} % Look at the grid
        %
        % Creates MGRID, 1-dim, 10 pts.
        %    >> mGrid('nPoints',10)
        %
        % Creates MGRID, 2-dim. In the object space:
        % 10 pts. in the first dimension linearly spaced between 0 and 20,
        % 45 pts in the second dimension linearly spaced between 5 and 40.
        % In the Fourier space: 
        % 10 pts. in the first dimension linearly spaced between 1 and 10,
        % 45 pts. in the first dimension linearly spaced between 1 and 45
        %    >> a = mGrid('nPoints',[10 45], ...
        %                'xMin',[0 5], 'xMax', [20 40], ...
        %                'kMin', [1 1], 'kMax', [10 45]) 
        %
        % Second method examples:
        % (The argument 'kAxes' is optional, automatically set to
        % span from 1 to obj.dimensions in each dimension
        %
        % Creates MGRID, 2-dim, with the grid {[1 3 5 7 9], [1 4 7 10]}
        %    >> a = mGrid('xAxes', {1:2:10, 1:3:10})
        %
        % Creates MGRID, 1-dim, with the grid {[1 3 5 7 9]} in the
        % object space and with the grid {[11 13 15 17 19]} in the
        % Fourier space. (Note that both grids have same length.)
        %    >> a = mGrid( 'xAxes', {1:2:10}, 'kAxes', {11:2:20})

            p = inputParser;

            % First method arguments
            default_dimension = 1;
            default_nPoints = 200;
            
            default_xMin = [0];
            default_xMax = [2.0 * pi];
            
            default_kMin = [1];
            default_kMax = [default_nPoints];

            addParameter(p,'nPoints',default_nPoints,@isnumeric)
            addParameter(p,'xMin',default_xMin,@isnumeric)
            addParameter(p,'xMax',default_xMax,@isnumeric)
            addParameter(p,'kMin',default_kMin,@isnumeric)
            addParameter(p,'kMax',default_kMax,@isnumeric)

            % Second method arguments
            default_xAxes = 0; % Dummy default values
            defoult_kAxes = 0; % Explicitely redefined later 
            addParameter(p,'xAxes',default_kMin)
            addParameter(p,'kAxes',default_kMax)

            parse(p,varargin{:});

            % Figure out which arguments were passed and which
            % construction method will be used
            str_varargin = varargin(1:2:size(varargin,2));
            if any(strmatch('xAxes',str_varargin)) | ...
                    any(strmatch('kAxes',str_varargin))
                construction_method = 'second_method';
            else
                construction_method = 'first_method';
            end
            
            if strcmp(construction_method, 'first_method')
                % Set the dimension, depending on which argument is
                % specified
                if nargin ~= 0
                    if any(strmatch('nPoints',str_varargin))
                        obj.dimension = size(p.Results.nPoints,2);
                    elseif any(strmatch('xMin',str_varargin))
                        obj.dimension = size(p.Results.xMin,2);
                    elseif any(strmatch('xMax',str_varargin))
                        obj.dimension = size(p.Results.xMax,2);
                    elseif any(strmatch('kMin',str_varargin))
                        obj.dimension = size(p.Results.kMin,2);
                    elseif any(strmatch('kMax',str_varargin))
                        obj.dimension = size(p.Results.kMax,2);
                        obj.kMax = p.Results.kMax;
                    end
                    
                    % If explicitely stated, set kMin and kMax
                    if any(strmatch('kMin',str_varargin)) & ...
                            any(strmatch('kMax',str_varargin))
                        obj.kMin = p.Results.kMin;                    
                        obj.kMax = p.Results.kMax;
                    else
                        error(['Both (or none) of the arguments kMin, ' ...
                               'kMax must be specified.'])
                    end
                else
                    obj.dimension = default_dimension;
                    obj.kMin = p.Results.kMin;
                    obj.kMax = p.Results.kMax;
                end
                
                % Set parsed arguments
                obj.nPoints = p.Results.nPoints;
                obj.xMin = p.Results.xMin;
                obj.xMax = p.Results.xMax;
                
                % Check that all input arguments have the same
                % dimension
                if ~(size(obj.nPoints)==size(obj.xMin) & ...
                     size(obj.xMin)==size(obj.xMax) & ...
                     size(obj.xMax)==size(obj.kMin) & ...
                     size(obj.kMin)==size(obj.kMax))
                    error(strcat('Input arguments passed to the ',... 
                                 ' mGrid constructor must have ',...
                                 ' identical shapes. Provided: ',...
                                 ' nPoints: ', mat2str(size(obj.nPoints)),...
                                 ' xMin: ', mat2str(size(obj.xMin)),...
                                 ' xMax: ', mat2str(size(obj.xMax)),...
                                 ' kMin: ', mat2str(size(obj.kMin)),...
                                 ' kMax: ', mat2str(size(obj.kMax)),...
                                 '.'))
                end
                
                % Set xAxes and kAxes
                obj.xAxes = {};
                obj.kAxes = {};
                for dim = 1 : 1 : obj.dimension
                    obj.xAxes = [obj.xAxes, linspace(obj.xMin(:,dim), ...
                                                     obj.xMax(:,dim), ...
                                                     obj.nPoints(:,dim))];
                    obj.kAxes = [obj.kAxes, linspace(obj.kMin(:,dim), ...
                                                     obj.kMax(:,dim), ...
                                                     obj.nPoints(:,dim))];
                end

            elseif strcmp(construction_method,'second_method')
                % Check that xAxes was provided
                if ~any(strmatch('xAxes',str_varargin))
                    error(['User input included optional parameter ' ...
                           'kAxes. This parameter makes input of ' ...
                           'xAxes necessary. Error: parameter xAxes' ...
                           'was not provided.'])
                end
                
                % Assure that xAxes is a horizontal vector
                tmp = size(p.Results.xAxes);
                if tmp(1) ~= 1
                    error(['xAxes must have dimension (1,n) for some ' ...
                           'natural number n. Provided dimensions:' ...
                           ''])
                end
                % Set xAxes
                obj.xAxes = p.Results.xAxes;
                tmp = size(obj.xAxes);
                obj.dimension = tmp(2);
                obj.nPoints = cellfun('length', p.Results.xAxes);
                
                obj.xMin = cellfun(@min, p.Results.xAxes);
                obj.xMax = cellfun(@max, p.Results.xAxes);
                
                % Check whether kAxes were provided; if not, create
                % defaults
                if any(strmatch('kAxes',str_varargin))
                    % Check that sizes of xAxes and kAxes match
                    if ~(all(cellfun('length', p.Results.xAxes) == ...
                             cellfun('length', p.Results.kAxes)))
                        error(['xAxes and kAxes size mismatch. ' ...
                               'Provided lengths: xAxes: ' ...
                               mat2str(cellfun('length', p.Results.xAxes))...
                               'kAxes: '...
                               mat2str(cellfun('length', p.Results.kAxes))...
                               '.'])
                    end
                    obj.kAxes = p.Results.kAxes;
                    obj.kMin = cellfun(@min, obj.kAxes);
                    obj.kMax = cellfun(@max, obj.kAxes);
                else
                    % See the end of the file for my_list function
                    obj.kAxes = myRange(cellfun('length', p.Results.xAxes));
                    obj.kMin = cellfun(@min, obj.kAxes);
                    obj.kMax = cellfun(@max, obj.kAxes);
                end
            
            end % End of the second method
        end % End of the constructor
    end % End of class methods
end % End of the class


% Small helpful function. Takes an array of natural numbers n1, n2, ...
% returns a cell array containing lists [1 2 3 ... n1], [1 2 3
% ... n2], ... etc.
function res = myRange(list_of_n)
    res = {};
    for  iList = 1 : 1 : length(list_of_n)
        res = [res, linspace(1, list_of_n(iList), list_of_n(iList))];
    end
end
