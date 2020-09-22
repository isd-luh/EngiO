% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Kursawe(x)
% KURSAWE is a two-objective, unconstrained, 3-dimensional test function.
%   [-5 -5 -5] <= x <= [5 5 5]
%
% LITERATURE:
%   F. Kursawe (1990) A variant of evolution strategies for vector
%   optimization,Proceedings of the International Conference on Parallel 
%   Problem Solving from Nature

if isempty(x)
    y=2;
    return
end

y= [0 0];

for i=1:2
    y(1)=y(1)-10*exp(-0.2*sqrt(x(i).^2+x(i+1).^2));
end

for i=1:3
    y(2)=y(2)+abs(x(i)).^0.8+5*sin(x(i).^3);
end

end