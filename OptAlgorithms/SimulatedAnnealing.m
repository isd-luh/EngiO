% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef SimulatedAnnealing < Optimizer
% SimulatedAnnealing is a single-objective Optimization Algorithm
% derived from Optimizer
%
% LITERATURE:
% S. Kirkpatrick (1983) Optimization by Simulated Annealing, Science
%
% K.-L. Du; M.N.S. Swamy (2016) Search and Optimization by Metaheuristics,
% Birkhäuser, Springer International Publishing AG Switzerland

    
    methods
        function this = SimulatedAnnealing()
            
            this.defaultparams = struct('CoreSize', 50, ...	% Number of observed cores in metal
                                        'T_0', 1, ...       % Initial value for metal temperature 
                                        'alpha', .95);      % Factor to generate new (lower) temperature 
        end
        
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            dim = numel(x_lb);
            
            % Initialize all states
            this.states = struct('x', zeros(this.params.CoreSize, dim),...           % Current core positions 
                                 'f', [],...                                         % Cache for function values 
                                 'Temp', this.params.T_0,...                         % Current temperature 
                                 'f_PersonalBest', Inf(this.params.CoreSize, 1),...  % Personal best function value 
                                 'f_GlobalBest', Inf,...                             % Global best function value 
                                 'numLoops', 0,...                                   % Current iteration step
                                 'dx', [],...                                        % Delta x for finding new neighbours
                                 'x_PersonalBest', [],...                            % Personal best core position
                                 'x_GlobalBest', zeros(1, dim));                     % Global best core position
                                
            % Initialize starting positions of cores (randomly)
            for i = 1:this.params.CoreSize
                if isempty(vartype)
                    this.states.x(i,:) = x_lb + (x_ub - x_lb) .* rand(size(this.states.x(i,:)));
                else
                    for j = 1:dim
                        if strcmpi(vartype{j}, 'real')
                            this.states.x(i, j) = x_lb(j) + (x_ub(j) - x_lb(j)) * rand;
                        elseif strcmpi(vartype{j}, 'int')
                            this.states.x(i, j) = randi([x_lb(j) x_ub(j)]);
                        end
                    end
                end
            end
            
            this.states.x_PersonalBest = this.states.x;
            
        end
        
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            dim = numel(x_lb);
            
            % Increment loop counter
            this.states.numLoops = this.states.numLoops+1;
        
            % First iteration: Go directly to terminate function
            if this.states.numLoops == 1
                samples = this.states.x;
                return 
            end      
        
            % All next iterations: Cool down temperature 
            this.states.Temp = this.states.Temp * this.params.alpha;

            % Find new neighbor for each core position 
            for i = 1:this.params.CoreSize
            
                this.states.dx(i,:) = (this.states.Temp/this.params.T_0) * ...
                    (x_lb +(x_ub - x_lb) .* rand(1, dim));

                % Calculate new states x + dx
                this.states.x(i,:) = this.states.x_GlobalBest + this.states.dx(i,:);
                
                % Limit positions to boundaries
                for j = 1:length(x_lb)
                    if this.states.x(i, j) < x_lb(j)
                        this.states.x(i, j) = x_lb(j);
                    elseif this.states.x(i, j) > x_ub(j)
                        this.states.x(i, j) = x_ub(j);
                    end
                end
                
                % Round to integer variables
                if ~isempty(vartype)
                    for j = 1:dim
                        if strcmpi(vartype{j}, 'int')
                            this.states.x(i, j) = round(this.states.x(i, j));
                        end
                    end
                end
                
            end
            
            % Store core positions as samples
            samples = this.states.x;
            
        end
        
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Add current results to cache for function values 
            this.states.f = objectiveValues;
            
            % Loop over all cores
            for j = 1:this.params.CoreSize
            
                % Determine acceptance of current results
                % If new minimum is found -> accept
                if objectiveValues(j) < this.states.f_PersonalBest(j)
                    this.states.f_PersonalBest(j) = objectiveValues(j);
                    this.states.x_PersonalBest(j,:) = samples(j,:);
                                
                % If current result > old result -> accept with certain probability
                else
                    DeltaEnergy = objectiveValues(j) - this.states.f_PersonalBest(j); % Difference between the two minimal values 
                    p = exp(-DeltaEnergy/this.states.Temp);                           % Probability factor
                
                    % If random value between 0 and 1 < p -> accept
                    if rand < p
                        this.states.f_PersonalBest(j) = objectiveValues(j);
                        this.states.x_PersonalBest(j,:) = samples(j,:);
                    end
                
                end
                
                % If best personal result is also best global: adjust
                if this.states.f_PersonalBest(j) < this.states.f_GlobalBest
                    this.states.f_GlobalBest = this.states.f_PersonalBest(j);
                    this.states.x_GlobalBest = this.states.x_PersonalBest(j,:);
                end

            end

            this.states.x = samples;
        
            terminate = false;
    
        end      
        
    end
end

