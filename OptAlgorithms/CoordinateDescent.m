% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef CoordinateDescent < Optimizer
    % COORDINATEDESCENT is a single-objective local Algorithm derived from Optimizer.
    %
    % LITERATURE:
    %   S.J. Wright (2015) Coordinate descent algorithms,
    %   Mathematical Programming
    
    methods
        
        function this = CoordinateDescent(states_init, options)
            
            this.defaultparams = struct('stepSize', 2.0);% Step size
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Initialize all states
            this.states = struct('nEvals', 0, ...
                'fBest', Inf, ...
                'sigma', zeros(length(x_lb), 1), ...
                'mean', zeros(length(x_lb), 1), ...
                'Idx', 0);
            
            N= length(x_lb); % Dimension
            
            XminPrime = zeros(N, 1);   % x^min
            XmaxPrime = zeros(N, 1);   % x^max
            xmean = zeros(N, 1);       % m
            sigma = zeros(N, 1);       % sigma
            
            % Loop over dimension
            for i = 1:N
                XminPrime(i, 1) = x_lb(i);
                XmaxPrime(i, 1) = x_ub(i);
                xmean(i, 1) = XminPrime(i, 1) + (rand())*(XmaxPrime(i, 1) - XminPrime(i, 1));
                sigma(i, 1) = (XmaxPrime(i, 1) - XminPrime(i, 1)) / 4;
            end
            
            this.states.sigma = sigma;
            this.states.mean = xmean;
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % B search direction
            B= eye(length(x_lb));
            
            % Choose index i_x
            this.states.Idx = mod(this.states.Idx, length(x_lb))+1; % Cyclic
            %             this.states.Idx = randi(length(x_lb)); % Random
            
            % Assign x' (= zeros beside i_x-th item)
            dx = this.states.sigma(this.states.Idx, 1)*B(:, this.states.Idx);
            
            % Sample two candidate solutions
            x1 = this.states.mean-dx;
            x2 = this.states.mean+dx;
            
            samples = [x1, x2]'; %(2 x dim)
            
            % If candidates violate the boundary than limit its value to
            % Boundaries
            if x1(this.states.Idx)<x_lb(this.states.Idx)
                x1(this.states.Idx) = x_lb(this.states.Idx);
            elseif x1(this.states.Idx)>x_ub(this.states.Idx)
                x1(this.states.Idx) = x_ub(this.states.Idx);
            end
            if x2(this.states.Idx)<x_lb(this.states.Idx)
                x2(this.states.Idx) = x_lb(this.states.Idx);
            elseif x2(this.states.Idx)>x_ub(this.states.Idx)
                x2(this.states.Idx) = x_ub(this.states.Idx);
            end
            
            samples = [x1, x2]'; % (2 x dim)
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Check if new best objective value
            lsucc = 0;
            if objectiveValues(1)<this.states.fBest
                this.states.fBest = objectiveValues(1);
                this.states.mean = samples(1,:)';
                lsucc = 1;
            end
            if objectiveValues(2)<this.states.fBest
                this.states.fBest = objectiveValues(2);
                this.states.mean = samples(2,:)';
                lsucc = 1;
            end
            
            % Adapt step-size sigma depending on the success/unsuccess of the previous search
            if lsucc==1 % Increase the step-size
                this.states.sigma(this.states.Idx, 1) = this.states.sigma(this.states.Idx, 1)*this.params.stepSize;
            else % Decrease the step-size
                this.states.sigma(this.states.Idx, 1) = this.states.sigma(this.states.Idx, 1)*(1./this.params.stepSize);
            end
            
            terminate = false;
            
        end
        
    end
    
end