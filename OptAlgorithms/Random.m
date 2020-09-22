% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef Random < Optimizer
    % RANDOM is a Monte-Carlo sampling implementation derived from Optimizer.
    %
    
    methods
        function this = Random(states_init, options)
            
            % default parameters
            this.defaultparams = struct('batchSize', 100); % Number of samples generated in one go
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % there are no states
            this.states = struct();
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % get number of design variables
            nDim = size(x_lb, 2);
            
            % Generate random samples
            vRandom = rand(this.params.batchSize, nDim);
            
            % scale to design space
            for iDim = 1:nDim
                samples(:, iDim) = vRandom(:, iDim) .* (x_ub(:, iDim) - x_lb(:, iDim)) + x_lb(:, iDim);
            end
            
            % Integer variables are rounded
            if ~isempty(vartype)
                for m = 1:nDim
                    if strcmpi(vartype{m}, 'int')
                        samples(:, m) = round(samples(:, m));
                    end
                end
            end
            
        end
        
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % termination by framework
            terminate = false;
            
        end
        
    end
end