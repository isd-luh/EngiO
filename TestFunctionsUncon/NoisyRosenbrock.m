% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = NoisyRosenbrock(x)
% NOISYROSENBROCK is a single-objective, unconstrained, i-dimensional test function.
%   Calls the ROSENBROCK function and adds noise.
%   Global minimum for i=... is at f(1,...,1)=0.
%   [-30 ... -30] <= x <= [30 ... 30]
%
%   See also ROSENBROCK.

if isempty(x)
    y=1;
    return
end

noise = 20*(prod(sin(20.0*x)) + 1.0);

y = Rosenbrock(x) + noise;
end