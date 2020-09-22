% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef Harmony < Optimizer
    % HARMONY is a single-objective Harmony Search algorithm derived from
    % Optimizer.
    %
    % LITERATURE:
    %   Z. W Geem, J. H. Kim, and G. V. Loganathan (2001) A new heuristic
    %   optimization algorithm: harmony search, Simulation 76.2: pp. 60-68
    
    methods
        
        function this = Harmony()
            
            this.defaultparams = struct('popSize', 25,... % Population size
                'HMCR', 0.90,...  % Harmony Memory Consideration Rate
                'PAR', 0.50,...   % Pitch Adjusting Rate
                'bw', 0.01);      % Bandwidth of change at PAR.
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Initialize all states
            this.states = struct('f',  ones (this.params.popSize, 1)*inf,...       % Objective function values for solutions inside HM
                'HM', zeros(this.params.popSize, length(x_lb)) ); % Harmony Memory
            
            % Initialize Harmony Memory
            for j = 1:this.params.popSize
                if isempty(vartype)
                    this.states.HM(j,:) = x_lb+(x_ub-x_lb).*rand(size(x_lb));
                else
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'real')
                            this.states.HM(j, m) = x_lb(m)+(x_ub(m)-x_lb(m)).*rand(1);
                        elseif strcmpi(vartype{m}, 'int')
                            this.states.HM(j, m) = randi([x_lb(m) x_ub(m)]);
                        end
                    end
                end
            end
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % initialize new sample
            xnew = zeros(1, length(x_lb));
            
            % assign sample coordinates
            for j = 1:size(xnew, 2)
                
                if rand <= this.params.HMCR % take val from inside HM
                    index = ceil(rand*this.params.popSize); % randomize index of target sol inside HM
                    val = this.states.HM(index, j); % get the val from HM
                    if rand <= this.params.PAR % check if a small shift is applicable
                        valBW = this.params.bw * (x_ub(j) - x_lb(j)); % calculate shift window
                        dval = (randn+0.5)*valBW - 0.50*valBW; % center BW around 0%
                        val = val +dval; % apply shift
                    end
                else % if not taken from inside HM
                    val = x_lb(j)+(x_ub(j)-x_lb(j))*rand; % just randomly assume it.
                end
                % Check boundaries
                val = min(max(x_lb(j), val), x_ub(j));
                xnew(j) = val;
                
                % Integer variables are rounded
                if ~isempty(vartype)
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'int')
                            xnew(m) = round(xnew(m));
                        end
                    end
                end
            end
            
            % Job done, return.
            samples = xnew;
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % find worst sample in buffer
            [maxVal, maxIdx] = max(this.states.f);
            
            % replace with current sample if better
            if(objectiveValues < maxVal)
                this.states.HM(maxIdx,:) = samples;
                this.states.f(maxIdx) = objectiveValues;
            end
            
            % termination by framework
            terminate = false;
        end
        
    end
    
end