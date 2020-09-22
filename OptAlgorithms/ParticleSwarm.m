% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef ParticleSwarm < Optimizer
    % PARTICLESWARM is a single-objective Particle Swarm Algorithm derived from
    % Optimizer.
    %
    % LITERATURE:
    %   J. Kennedy and R. Eberhart (1995) Particle Swarm
    %   Optimization, Proceedings of the IEEE international conference on neural
    %   networks
    
    methods
        
        function this = ParticleSwarm()
            
            this.defaultparams = struct('c1', 2.0,...       % First acceleration factor that considers the impact of personal best solution
                'c2', 2.0,...                             % Second acceleration factor that considers the impact of global best solution
                'c3', 0.0,...                             % Third acceleration factor that considers the impact of passive congregation
                'omega1', 0.9,...                         % Initial inertia weight
                'omega2', 0.4,...                         % Final inertia weight (linear interpolation between values)
                'popSize', 50);                           % Population size
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Initialize all states
            this.states = struct('f', zeros(this.params.popSize, 1),... % Objective function values for current particle positions
                'f_p', inf(this.params.popSize, 1),... % Objective function values for personal best particle positions
                'f_g', inf,... % Objective function values for global best particle position
                'numLoops', 0,... % Current generation
                'v', zeros(this.params.popSize, length(x_lb)),... % Particle velocities
                'w', this.params.omega1*ones(1, this.params.popSize),... % Inertia weight factors
                'x', zeros(this.params.popSize, length(x_lb)),... % Current particle positions
                'x_p', zeros(this.params.popSize, length(x_lb)),... % Personal best particle positions
                'x_g', zeros(1, length(x_lb))); % Global best particle positions
            
            this.states.loopOmegaMax = min(this.options.maxIters, this.options.maxEvals/this.params.popSize); % Iteration loop where omega value reaches maximumloopOmegaMax = min(this.options.maxLoops, this.options.maxEvals/popSize); % Iteration loop where omega value reaches maximum
            
            % Initialize particle positions
            for j = 1:this.params.popSize
                if isempty(vartype)
                    this.states.x(j,:) = x_lb+(x_ub-x_lb).*rand(size(this.states.x(j,:)));
                else
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'real')
                            this.states.x(j, m) = x_lb(m)+(x_ub(m)-x_lb(m)).*rand(1);
                        elseif strcmpi(vartype{m}, 'int')
                            this.states.x(j, m) = randi([x_lb(m) x_ub(m)]);
                        end
                    end
                end
            end
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % Increment loop counter
            this.states.numLoops = this.states.numLoops+1;
            if this.states.numLoops==1
                samples = this.states.x;
                return
            end
            
            % Loop over all particles
            for j = 1:this.params.popSize
                
                % Calculate inertia weight (approach by Shi and Eberhart)
                this.states.w = this.params.omega1+(this.params.omega2-this.params.omega1)/(this.states.loopOmegaMax-1)*(this.states.numLoops-1);
                
                % Calculate new particle velocities
                this.states.v(j,:) = this.states.w*this.states.v(j,:)... % Inertia term
                    +this.params.c1*rand()*(this.states.x_p(j,:)-this.states.x(j,:))... % Personal compulsion term
                    +this.params.c2*rand()*(this.states.x_g-this.states.x(j,:))... % Global compulsion term
                    +this.params.c3*rand()*(this.states.x(randi(this.params.popSize))-this.states.x(j,:)); % Passive congregation term % Passive congregation term
                
                % See: 10.1016/j.biosystems.2004.08.003
                % We propose a hybrid PSO with passive congregation: where Ri
                % is a particle selected randomly from the swarm, c3 the
                % passive congregation coefficient ...
                
                % Integer variables are rounded
                if ~isempty(vartype)
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'int')
                            this.states.v(j, m) = round(this.states.v(j, m));
                        end
                    end
                end
                
                this.states.x(j,:) = this.states.x(j,:)+this.states.v(j,:);
                
                % If one particle violates the boundary than limit its position
                % and modify the velocity by reflecting
                for m = 1:length(x_lb)
                    if this.states.x(j, m)<x_lb(m)
                        this.states.x(j, m) = x_lb(m);
                        this.states.v(j, m) =-0.5*this.states.v(j, m);
                    elseif this.states.x(j, m)>x_ub(m)
                        this.states.x(j, m) = x_ub(m);
                        this.states.v(j, m) =-0.5*this.states.v(j, m);
                    end
                end
                
            end
            
            % Store particle positions as samples
            samples = this.states.x;
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Loop over all particles
            for j = 1:this.params.popSize
                
                % Check for all particles whether current solution is best personal solution
                if objectiveValues(j)<this.states.f_p(j)
                    this.states.f_p(j) = objectiveValues(j);
                    this.states.x_p(j,:) = samples(j,:);
                    % Check whether new best personal solution is also best global solution
                    if this.states.f_p(j)<this.states.f_g
                        this.states.x_g = this.states.x_p(j,:);
                        this.states.f_g = this.states.f_p(j);
                    end
                end
                
            end
            
            this.states.f = objectiveValues;
            
            terminate = false;
            
        end
        
    end
    
end