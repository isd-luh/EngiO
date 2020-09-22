% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function res = TenBarTruss(x)
% TENBARTRUSS is a single-objective, constrained, 10-dimensional test function.
%   The design variables specify the cross-sectional areas of the ten
%   members.
%   [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1] <= x <= 
%   [33.5 33.5 33.5 33.5 33.5 33.5 33.5 33.5 33.5 33.5]
%   
%   There exist different, possible sets of constraints for the TENBARTRUSS
%   test function. In this example, only the stress is constrained.
%
% LITERATURE:
%   Zabinsky, Z. B. (2003) Stochastic Adaptive Search for Global
%   Optimization, Springer Science+Business Media, LLC

if isempty(x)
    res={1, 1, 0};
    return
end

% -------------------------------------------------------------------------
% (1) Geometry
% -------------------------------------------------------------------------

sqin2sqm = 0.0254^2; % inch^2 to m^2
A = x*sqin2sqm;

% Node coordinates
len = 9.144;         % (m)
node = [2*len len    % node 1
        2*len 0      % node 2
        len   len    % node 3
        len   0      % node 4
        0     len    % node 5
        0     0];    % node 6
    
% Number of nodes
num_nodes = size(node,1);

% Elements
elem = [3 5          % elem 1
        1 3          % elem 2
        4 6          % elem 3
        2 4          % elem 4
        3 4          % elem 5
        1 2          % elem 6
        4 5          % elem 7
        3 6          % elem 8
        2 3          % elem 9
        1 4];        % elem 10
    
% Number of elements
num_ele = size(elem,1);    

% Global degrees of freedom allocated to elements  
ele_dof = [5 6 9 10     % elem 1
           1 2 5 6      % elem 2
           7 8 11 12    % elem 3
           3 4 7 8      % elem 4
           5 6 7 8      % elem 5
           1 2 3 4      % elem 6
           7 8 9 10     % elem 7
           5 6 11 12    % elem 8
           3 4 5 6      % elem 9
           1 2 7 8];    % elem 10
       
% Number of degrees of freedom
num_dof = max(max(ele_dof));

% -------------------------------------------------------------------------
% (2) Material
% -------------------------------------------------------------------------

% Young's modulus
E(1:10) = 6.89476*10^10; %(Pa / N/m^2)

% -------------------------------------------------------------------------
% (3) Loads & boundary conditions
% -------------------------------------------------------------------------

% Applied load at degrees of freedom 
force = zeros(2*num_nodes,1);
force(4) = -444822.16;      % (N)
force(8) = -444822.16;      % (N) 

% Boundary conditions at degrees of freedom 
dofBOUND = [9, 10, 11, 12];
displacement = zeros(num_ele,1);
for e = dofBOUND(1):dofBOUND(end)
    displacement(e,1) = 0;
end

% -------------------------------------------------------------------------
% (4) Calculation of stiffness matrix 
% -------------------------------------------------------------------------

stiffness = zeros(num_dof);
L = zeros(num_ele,1);
C = zeros(num_ele,1);
S = zeros(num_ele,1);
for e = 1:num_ele
    L(e) = sqrt((node(elem(e,2),1)-node(elem(e,1),1))^2+(node(elem(e,2),2)-node(elem(e,1),2))^2);
    C(e) = (node(elem(e,2),1)-node(elem(e,1),1))/L(e);
    S(e) = (node(elem(e,2),2)-node(elem(e,1),2))/L(e);
    k = (A(e)*E(e)/L(e)*[C(e)*C(e) C(e)*S(e) -C(e)*C(e) -C(e)*S(e);C(e)*S(e) S(e)*S(e) -C(e)*S(e) -S(e)*S(e);...
        -C(e)*C(e) -C(e)*S(e) C(e)*C(e) C(e)*S(e); -C(e)*S(e) -S(e)*S(e) C(e)*S(e) S(e)*S(e)]);
        
    % Extract the rows of ele_dof for each element e
    ele_dof_vec = ele_dof(e,:);
    for i = 1:4
        for j = 1:4
            stiffness(ele_dof_vec(i),ele_dof_vec(j)) = stiffness(ele_dof_vec(i),ele_dof_vec(j))+k(i,j);
        end
    end
end

% -------------------------------------------------------------------------
% (5) Calculation of displacement vector
% -------------------------------------------------------------------------

num_dofUNBOUND = 8;
displacement(1:num_dofUNBOUND) = stiffness(1:num_dofUNBOUND,1:num_dofUNBOUND)\force(1:num_dofUNBOUND);

% Calculation force for boundary
for i = 1:size(dofBOUND,1)
    force(dofBOUND(i),1) = stiffness(dofBOUND(i),:)*displacement;
end

% -------------------------------------------------------------------------
% (6) Calculation of stress in all elements
% -------------------------------------------------------------------------

sigma = zeros(num_ele,1);
for e = 1:num_ele
    sigma(e) = (E(e)/L(e))*[-C(e) -S(e) C(e) S(e)]*displacement((ele_dof(e,:))');
end

% -------------------------------------------------------------------------
% (7) Inequality constraint for stress 
% -------------------------------------------------------------------------
 
sigma_lim = 1.72369*10^8*ones(num_ele,1); % (Pa / N/m^2)

% Possibility to further constrain stress of truss member 9
% sigma_lim(9) = 5.17107*10^8; 

% Constraint: no member can be overused
diff_sigma = sigma_lim - abs(sigma);         
g = min(diff_sigma);                      % if g < 0: not satisfyed

% -------------------------------------------------------------------------
% (7) Weight of the structure
% -------------------------------------------------------------------------

dens = 2767.99;           % (kg/m^3)
elem_weight = L.*A'*dens; % (kg)
y = sum(elem_weight);     % (kg)

% -------------------------------------------------------------------------
res={y, g}; 
end