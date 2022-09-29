[input_pos, ~, ~] = xlsread('input_pos.xlsx');
[input_load, ~, ~] = xlsread('input_load.xlsx');
[input_hinge, ~, ~] = xlsread('input_hinge.xlsx');

%% columns of input_pos corresponds to
% element number =1
% i = 2 
% j = 3
% xi = 4
% yi = 5
% xj = 6
% yj = 7
%% columns of input_load corresponds to
% element number =1
% Fx = 2 
% Fy = 3
%% columns of input_hinge corresponds to
% element number =1
% Fx = 2 
% Fy = 3
%% Programme
le = sqrt((input_pos(:,6)-input_pos(:,4)).^2+(input_pos(:,7)-input_pos(:,5)).^2);
cos = (input_pos(:,6)-input_pos(:,4))./le;
sin = (input_pos(:,7)-input_pos(:,5))./le;
nel = length(input_pos);
nnodes = length(input_hinge);
A = 2.5 * 10^(-3 );
E = 2* 10^7; 

%local stiffness matrix
k = zeros(4,4,nel);
%transformation matrix
T = zeros(4,4,nel);
%Lather matrix
L = zeros(4,nnodes*2,nel);
%connectivity matrix
C = [];

for e = 1:nel
    
    k(:,:,e) = (A*E/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
    T(:,:,e) = [cos(e) sin(e) 0 0 ; -sin(e) cos(e) 0 0 ; 0 0 cos(e) sin(e) ; 0 0 -sin(e) cos(e)];
    i = input_pos(e,2);
    j = input_pos(e,3);
    L(1,((2*i)-1),e) = 1;
    L(2,(2*i),e) = 1;
    L(3,((2*j)-1),e) = 1;
    L(4,(2*j),e) = 1;
    C(e,1) = 2*i-1;
    C(e,2) = 2*i;
    C(e,3) = 2*j-1;
    C(e,4) = 2*j;
end

%% Assembly
K = zeros(2*nnodes,2*nnodes);
for e = 1:nel
    K = K + ((L(:,:,e))'*(T(:,:,e))'*k(:,:,e)*T(:,:,e)*L(:,:,e));
end

% Force vector
F = zeros(2*nnodes,1);
for i = 1:nnodes
    F(2*i-1) = input_load(i,2);
    F(2*i) = input_load(i,3);
end

% Penalty Approach
K_hat = K;
F_hat = F;
k_mean = 0;
for i = 1:20
    k_mean = k_mean + K(i,i); 
end
beta = k_mean * (10^7);
%hinge at position 1,2,12
for i = 1:nnodes
    if input_hinge(i,2)== 1
        K_hat((2*i-1),(2*i-1)) = beta + K(2*i-1);
        F_hat(2*i-1) = beta*0;
%         K_hat((2*i-1),(2*i-1)) = 1;
%         F_hat = F_hat - (0.*K((2*i-1),:))';
    end
    if input_hinge(i,3)== 1
        K_hat ((2*i),(2*i)) = beta + K(2*i);
        F_hat(2*i) = beta*0;
%         K_hat((2*i),(2*i)) = 1;
%         F_hat =  F_hat - (0.*K((2*i),:))';
    end
end
% Solution  
global_dof = K_hat\F_hat;
for i = 1:nnodes
    if input_hinge(i,2)== 1
        global_dof((2*i-1)) = 0;
    end
    if input_hinge(i,3)== 1
        global_dof((2*i)) = 0;
    end
end
%% Post Processing
local_dof = zeros(4,nel);
ax_force = zeros(nel,1);
for e = 1:nel
    local_dof(:,e) = T(:,:,e)*L(:,:,e)*global_dof;
    ax_force (e,1)= ((local_dof(3,e)-local_dof(1,e))*E*A)/(le(e));
end
    
    


Disp = cell(nnodes+1,3);
Disp(1,:) = {'Node_number','Displacement_x (in m)','Displacement_y (in m)'};
for i = 1:nnodes
    Disp(i+1,1) = num2cell(i);
    Disp(i+1,2) = num2cell(global_dof((2*i)-1));
    Disp(i+1,3) = num2cell(global_dof((2*i)));
end
Ax_force = cell(nel+1,2);
Ax_force(1,:) = {'Element_Number','Axial_Force (in N)'};
for i = 1:nel
    Ax_force(i+1,1) = num2cell(i);
    Ax_force(i+1,2) = num2cell(ax_force(i));
end

xlswrite('output_disp.xlsx', Disp);
xlswrite("output_ax_force.xlsx", Ax_force);
%% Plotting
for e = 1:nel
    plot([input_pos(e,4),input_pos(e,6)],[input_pos(e,5),input_pos(e,7)],'k');
    axis equal;
    hold on;
end
for e = 1:nel
    i = input_pos(e,2);
    j = input_pos(e,3);
    plot([(input_pos(e,4)+global_dof(2*i-1)),(input_pos(e,6)+global_dof(2*j-1))],[(input_pos(e,5)+global_dof(2*i)),(input_pos(e,7)+global_dof(2*j))],'r');
    hold on;
end



