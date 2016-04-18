classdef phaseRetrievalAlgorithm
%  phaseRetrievalAlgorithm - class to store appoximations and error data
%
% Class attributes:
%   obj.initial  - instance of molecule, stores initial
%                  approximation
%   obj.current  - instance of molecule, stores current
%                  approximation
%   obj.solution - instance of molecule, stores ideal result
%                  (if known)
%   obj.energy   - array storing energies corresponding to all
%                  approximations
%   obj.num_steps -total number of already performed steps
%
% Class methods:
%   obj.update   - perform an update step in the algorithm. 
%
% For examples, type 'doc phaseRetrievalAlgorithm' and see
% Constructor info.
% <UNFINISHED>
    properties
        initial % Initial approximation
        current % Current approximation
        update_params % Parameters for the update routine
    end
    properties (SetAccess=private)
        ftmod     % Known Fourier modulus (Measured sqrt of
                  % intensity)

        energy    % Energy vector
        num_steps % Total number of performed steps

        solution  % Exact solution (if known)

        update_function;    % Handle to the update function
    end
    methods
        function obj = phaseRetrievalAlgorithm(update_function_, ...
                                               solution_, update_params_)
            % Class constructor
            %
            % solution_ may be either an instance of the class
            % molecule, or the Fourier modulus (square root of the
            % measured signal.
            obj.update_function = update_function_;
            
            if isa(solution_, 'molecule')
                obj.solution = solution_;
                obj.ftmod = abs(obj.solution.ft);

                % Set the initial approximation
                disp('Info: Using initial phase 0. You may set phase using')
                disp('>> obj = obj.initial.set_phase(...)');
                obj.initial = molecule('ft', obj.ftmod, ...
                                       'grid', obj.solution.grid);
               
                disp(['Info: the exact solution of the problem was' ...
                      ' specified.']);
            else
                obj.solution = 'None';
                obj.ftmod = solution_;

                % Set the initial approximation
                disp('Info: Using initial phase 0. You may set phase using')
                disp('>> obj = obj.initial.set_phase(...)');
                obj.initial = molecule('ft', obj.ftmod);

                disp(['Info: the exact solution of the problem was' ...
                      ' not specified. Proceeding with FT modulus...']);
            end
            
            obj.current = obj.initial;
            obj.energy = [];
            obj.num_steps = 0;
            
            if nargin == 2
                obj.update_params = 'None';
            else
                obj.update_params = update_params_;
            end
        end
        
        function obj = update(obj, num_updates)
            if obj.update_params == 'None'
                for i = 1 : 1 : num_updates
                    [new_density, new_energy] = ...
                        obj.update_function(obj.current.density, obj.ftmod);
                    obj.current = obj.current.set_density(new_density);
                    obj.energy = [obj.energy new_energy];
                end
            else
                for i = 1 : 1 : num_updates
                    [new_density, energy] = ...
                        obj.update_function(obj.current.density, ...
                                            obj.ftmod, update_params);
                    obj.current = obj.current.set_density(new_density);
                    obj.energy = [obj.energy new_energy];
                end
            end
        end % update(obj, num_updates)

        % Longer functions are specified in their own files
        plot(obj, varargin);
    end
end