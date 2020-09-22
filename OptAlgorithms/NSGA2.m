% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef NSGA2 < Optimizer
    % NSGA2 is a multi-objective Genetic Algorithm derived from Optimizer.
    %
    % LITERATURE:
    %   K. Deb, S. Agrawal, A. Pratapj, T. Meyarivan (2000) A fast elitist
    %   non-dominated sorting genetic algorithm for multi-objective
    %   optimization: NSGA-II, International conference on parallel problem
    %   solving from nature, Springer
    methods
        
        function this = NSGA2()
            
            this.defaultparams = struct('eta_c', 10, ...% Distribution index for real variable SBX crossover
                'eta_m', 10, ...                        % Distribution index for real variable polynomial mutation
                'popSize', 500, ...                     % Population size
                'pcross_real', 0.9,...                  % Probability of crossover of real variable
                'pmut_real', 0.33...                    % Probability of mutation of real variable
                );
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Check for popsize
            if ( mod(this.params.popSize, 4) ~= 0 )
                error('Population of %d is not divisible by 4!', this.params.popSize);
            end
            
            % Initialize all states
            this.states = struct('f', inf(this.params.popSize, nObj),...    % Objective function values for current samples
                'f_g', inf(this.params.popSize, nObj),...                   % Objective function values for global best sample
                'f_lastGen', inf(this.params.popSize, nObj),...             % Current generation
                'x', zeros(this.params.popSize, length(x_lb)),...           % Current samples
                'x_lastGen', zeros(this.params.popSize, length(x_lb)), ...  % Samples from last iteration loop
                'initialized', false);                                      % Flag for initialization
            
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
            if ~this.states.initialized
                samples = this.states.x;
                return
            end
            
            % Allocate states
            this.states.x_lastGen = this.states.x; % copy processed samples
            this.states.f_lastGen = this.states.f;
            
            % Non-dominated sorting
            [cFront, vCountFront]= NonDomSort(this);
            % Sort by crowding distance and get crowding distances
            Crowd_Dist = cell([1 length(vCountFront)]);
            for iFront = 1:length(vCountFront)
                [cFront, Crowd_Dist]= CrowDistSort(this, nObj, cFront, vCountFront, iFront, Crowd_Dist);
            end
            
            % Get crowding distances per individual
            vCrowdDists = nan(this.params.popSize, 1);
            for iFront = 1:length(vCountFront)
                vCrowdDists(cFront{iFront}) = Crowd_Dist{iFront};
            end
            
            % Matrix for new population
            new_pop = zeros(this.params.popSize, numel(x_lb));
            
            % Random permutation for selection
            a1 = randperm(this.params.popSize);
            a2 = randperm(this.params.popSize);
            
            for j = 1:4:this.params.popSize
                % Tournament selection of parents
                parent1 = this.tournament(a1(j)  , a1(j+1), vCrowdDists);
                parent2 = this.tournament(a1(j+2), a1(j+3), vCrowdDists);
                % Crossover
                [new_pop(j,:), new_pop(j+1,:)]= this.crossover(parent1, parent2, x_lb, x_ub); %parent = index of parent
                
                % Tournament selection of parents
                parent1 = this.tournament(a2(j)  , a2(j+1), vCrowdDists);
                parent2 = this.tournament(a2(j+2), a2(j+3), vCrowdDists);
                % Crossover
                [new_pop(j+2,:), new_pop(j+3,:)]= this.crossover(parent1, parent2, x_lb, x_ub);
            end
            
            % Samples is the new generation after mutation
            samples = this.mutation(new_pop, x_lb, x_ub);
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            nObj = size(objectiveValues, 2);
            terminate = false;
            
            if ~this.states.initialized
                this.states.f = objectiveValues;
                this.states.x = samples;
                
                this.states.initialized = true;
                return
            else
                % Combine parent and children population (R_t = P_t union Q_t)
                this.states.f = [this.states.f_lastGen;objectiveValues];
                this.states.x = [this.states.x_lastGen;samples];
            end
            
            % Non-dominated sorting
            [cFront, vCountFront]= NonDomSort(this);
            
            % Choose fronts for parent population
            vCumCountFront = cumsum(vCountFront);
            for i = 1:numel(vCumCountFront)
                if vCumCountFront(i)>= this.params.popSize
                    iLastFront = i;
                    break
                end
            end
            
            % Crowding distance sorting
            for iFront = 1:iLastFront
                [cFront]= CrowDistSort(this, nObj, cFront, vCountFront, iFront);
            end
            
            % Fill parent population
            idxParents = NaN(this.params.popSize, 1);
            for i = 1:iLastFront
                if i==1
                    idxParents(1:vCumCountFront(i)) = cFront{i};
                elseif i == iLastFront
                    remain = this.params.popSize-vCumCountFront(i-1);
                    idxParents(vCumCountFront(i-1)+1:end) = cFront{i}(1:remain);
                else
                    idxParents(vCumCountFront(i-1)+1:vCumCountFront(i)) = cFront{i};
                end
            end
            
            this.states.x = this.states.x(idxParents,:);
            this.states.f = this.states.f(idxParents,:);
        end
        
    end
    
    % Private helper functions
    methods(Access = private)
        
        % Tournament selection
        function ind = tournament(this, ind1, ind2, vCrowdDists)
            % Check dominance
            if(all(this.states.f(ind1,:)<this.states.f(ind2,:)))
                ind = ind1;
                return
            end
            if(all(this.states.f(ind1,:)>this.states.f(ind2,:)) )
                ind = ind2;
                return
            end
            
            % Compare by crowding distance
            if vCrowdDists(ind1) > vCrowdDists(ind2)
                ind = ind1;
                return
            end
            if vCrowdDists(ind1) < vCrowdDists(ind2)
                ind = ind2;
                return
            end
            
            % Decide by chance
            if rand <= 0.5
                ind = ind1;
            else
                ind = ind2;
            end
        end
        
        
        % Crossover
        function [child1, child2]= crossover(this, parent1, parent2, x_lb, x_ub)
            child1 = this.states.x(parent1,:);
            child2 = this.states.x(parent2,:);
            
            % Check whether crossover will happen
            if rand > this.params.pcross_real
                return
            end
            
            % For each design variable
            for i = 1:numel(x_lb)
                % Crossover will not happen to every design variable
                if rand <=0.5
                    if this.states.x(parent1, i) ~= this.states.x(parent2, i)
                        % Sort parent coordinates
                        y1 = min(this.states.x(parent1, i), this.states.x(parent2, i) );
                        y2 = max(this.states.x(parent1, i), this.states.x(parent2, i) );
                        
                        % Crossover parameter
                        rnd = rand;
                        
                        beta = 1.0 + (2.0*(y1-x_lb(i))/(y2-y1));
                        alpha = 2.0 - (beta^(-(this.params.eta_c+1.0)));
                        if (rnd <= (1.0/alpha))
                            betaq = (rnd*alpha)^(1.0/(this.params.eta_c+1.0));
                        else
                            betaq = (1.0/(2.0 - rnd*alpha))^(1.0/(this.params.eta_c+1.0));
                        end
                        
                        % First crossover result
                        c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                        
                        beta = 1.0 + (2.0*(x_ub(i)-y2)/(y2-y1));
                        alpha = 2.0 - (beta^(-(this.params.eta_c+1.0)));
                        if (rnd <= (1.0/alpha))
                            betaq = (rnd*alpha)^(1.0/(this.params.eta_c+1.0));
                        else
                            betaq = (1.0/(2.0 - rnd*alpha))^(1.0/(this.params.eta_c+1.0));
                        end
                        
                        % Second crossover result
                        c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                        
                        % Clamp to design variable domain
                        c1 = max(min(c1, x_ub(i)), x_lb(i));
                        c2 = max(min(c2, x_ub(i)), x_lb(i));
                        
                        % Assign crossover result randomly to children
                        if (rand<=0.5)
                            child1(i) = c2;
                            child2(i) = c1;
                        else
                            child1(i) = c1;
                            child2(i) = c2;
                        end
                    end
                end
            end
        end
        
        % Mutation
        function [mutated_pop]= mutation(this, new_pop, x_lb, x_ub)
            
            mutated_pop = new_pop;
            
            for i = 1:this.params.popSize
                
                for j = 1:numel(x_lb)
                    if rand <= this.params.pmut_real
                        y = new_pop(i, j);
                        delta1 =(y-x_lb(j))/(x_ub(j)-x_lb(j));
                        delta2 =(x_ub(j)-y)/(x_ub(j)-x_lb(j));
                        rnd = rand;
                        mut_pow = 1.0/(this.params.eta_m+1.0);
                        if rnd<=0.5
                            xy  = 1.0-delta1;
                            val = 2.0*rnd+(1.0-2.0*rnd)*(xy^(this.params.eta_m+1.0));
                            deltaq = (val^mut_pow) - 1.0;
                        else
                            xy  = 1.0-delta2;
                            val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*((xy^(this.params.eta_m+1.0)));
                            deltaq = 1.0 - ((val^mut_pow));
                        end
                        
                        % Mutate variable
                        y = y + deltaq*(x_ub(j)-x_lb(j));
                        
                        % Clamp to design variable domain
                        y = max(min(y, x_ub(j)), x_lb(j));
                        
                        mutated_pop(i, j) = y;
                    end
                end
            end
            
        end
        
        function  [cFront, vCountFront]= NonDomSort(this)
            
            % Fast non-dominated sorting
            vCountFront = [];
            f = this.states.f;
            ind_org = [1:size(f, 1)]; % original indices
            ind_red = ind_org;       % reduced set of original indices
            i = 1;
            
            while sum(vCountFront)<size(f, 1)
                
                ind_front = this.paretoSort(f(ind_red,:)); % indices of front related to reduced set
                cFront{i}= ind_red(ind_front);
                vCountFront(i) = numel(cFront{i});
                ind_red = setdiff(ind_red, cFront{i});
                
                i = i+1;
            end
        end
        
        % Crowding Distance Sorting
        function  [cFront, Crowd_Dist]= CrowDistSort(this, nObj, cFront, vCountFront, iFront, Crowd_Dist)
            
            nFront = vCountFront(iFront); % number of solutions in front set
            idxComb = cFront{iFront};     % index of solutions in combined set R_t
            CrowdDist = zeros(nFront, 1);  % initialize crowding distance
            
            for m = 1:nObj
                [~, idxFront]= sortrows(this.states.f(idxComb, 1:nObj), m);     % sort using objective value
                CrowdDist(idxFront(1)) = inf;                        % so that boundary points are always selected
                CrowdDist(idxFront(nFront)) = inf;
                
                if numel(idxFront) <= 2
                    continue
                end
                
                j = 2:(nFront-1);
                CrowdDist((idxFront(j)))... % for all other points
                    = CrowdDist((idxFront(j)))...
                    +(this.states.f(idxComb(idxFront(j+1)), m)...
                    -this.states.f(idxComb(idxFront(j-1)), m));
            end
            
            [~, index] = sort(CrowdDist);
            
            % Sort crowding distances
            Crowd_Dist{iFront}= num2cell(CrowdDist);
            Crowd_Dist{iFront}= CrowdDist(flipud(index));
            % Sort indices
            cFront{iFront}= idxComb(flipud(index));
            
        end
        
    end
    
end


