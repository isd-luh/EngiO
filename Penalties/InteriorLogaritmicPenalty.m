% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef InteriorLogaritmicPenalty
    % INTERIORLOGARITMICPENALTY is a static constraint handling technique.
    % The penalty function has logaritmic slope in the feasible region.
    %
    % LITERATURE:
    %   C.A. Coello Coello (2002) Theoretical and numerical
    %   constraint-handling techniques used with evolutionary algorithms: a
    %   survey of the state of the art, Comput. methods Appl. Mech. Engrg.
    
    properties
        r;
    end
    
    methods
        
        function this = InteriorLogaritmicPenalty(r)
            
            % Call base class constructor
            this.r = r;
            
        end
        
        function penalty = penaltyFunction(this, nObj, g)
            
            % Return logarithmic penalty if feasible(g>=0), Inf if otherwise
            penalty = -this.r*log(g)';
            if g<0
                penalty = Inf(1, nObj);
            end
            
        end
        
    end
    
end