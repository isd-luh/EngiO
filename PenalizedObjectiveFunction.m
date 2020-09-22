% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef PenalizedObjectiveFunction < handle
    % PENALIZEDOBJECTIVEFUNCTION penalizes the objective function according to the
    % constraint handling techniques selected. The optimization problem
    % should be stated as follows:
    %   minimize f(x) 
    %   subject to  g(x) >=0, 
    %               h(x) =0, 
    %   for x_lb <= x <= x_ub
    
    properties
        
        fun; % Objective function
        nObj=0; % Number of Objectives
        
        eq; % Equality constraint penalizers
        ineq; % Inequality constraint penalizers
        
    end
    
    methods(Access=public)
        function this = PenalizedObjectiveFunction(problem, cInequalPenalizers, cEqualPenalizers)
            % PENALIZEDOBJECTIVEFUNCTION is the constructor.
            %
            % - problem (function handle)
            %   The constrained problem function to minimize
            % - cInequalPenalties (cell array of penalizers)
            %   The Penalizers for inequality constraints
            % - cEqualPenalties (cell array of penalizers)
            %   The Penalizers for equality constraints
            
            % Store objective function and penalizers
            this.fun = problem;
            this.ineq = cInequalPenalizers;
            this.eq   = cEqualPenalizers;
            
            % Initial call -> returns number of objectives and constraints
            if nargin(this.fun)==1
                result = this.fun([]);
            else
                result = this.fun([], 1);
            end
            
            % Get final number of objectives by calling child class
            this.nObj = result{1};
            
            nIneq = result{2};
            nEq = result{3};
            
            % Test if the number of constraints matches the number of penalizers
            if(numel(this.ineq) ~= nIneq || numel(this.eq) ~= nEq)
                error('Constraint numbers must be consistent');
            end

        end
        
        function y=evaluate(this, x, index)
            % EVALUATE the constraint equation and add penalty to objective.
            %
            % - x (real vector)
            %   The design variable vector
            % - index (integer)
            %   The index for parallelized evaluation
            %
            % OUTPUT:
            % - y (double scalar|vector)
            %   Penalized objective values
            
            % Handles initial call
            if isempty(x)
                y = this.nObj;
                return
            end
            
            % Evaluate objective function
            if nargin(this.fun)==1
                result = this.fun(x);
            else
                result = this.fun(x, index);
            end
            
            % Copy result to objectives
            f_out = result{1};
            
            % Get penalty vector by calling penalizer

            % Inequality
            for iIneq = 1:numel(this.ineq)
                pen = this.ineq{iIneq}.penaltyFunction(this.nObj, result{2}(iIneq));
                f_out = f_out + pen;
            end
            
            % Equality
            for iEq = 1:numel(this.eq)
                pen = this.eq{iEq}.penaltyFunction(this.nObj, -abs(result{3}(iEq)) );
                f_out = f_out + pen;
            end
                        
            % Return objective value with added penalty (expanded objective)
            y = f_out;
        end
        
    end
    
end