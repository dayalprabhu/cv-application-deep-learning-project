[input_pos, ~, ~] = xlsread('input_pos.xlsx');
[input_load, ~, ~] = xlsread('input_load.xlsx');
[input_hinge, ~, ~] = xlsread('input_hinge.xlsx');
[input_disp_b, ~, ~] = xlsread('input_disp_b.xlsx');

%% columns of input_pos corresponds to
% element number =1
% i = 2 
% j = 3
% xi = 4
% yi = 5
% xj = 6
% yj = 7
% type = 8
% if type == 1, element is bottom chord.
% if type == 2, element is top chord.
% if type == 3, element is vertical.
% if type == 4, element is diagonal.
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
unknown = 4;
%local stiffness matrix
k = zeros(4,4,nel);
%transformation matrix
T = zeros(4,4,nel);
%Lather matrix
L = zeros(4,nnodes*2,nel);
%connectivity matrix
C = [];

for e = 1:nel
    k(:,:,e) = (A/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
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

k_temp = zeros(4,4,nel);
for i = 1:nel
    k_temp(:,:,i) = T(:,:,i)'*k(:,:,i)*T(:,:,i);
end
res_DoFs = [];
for i = 1:length(input_hinge)
     if input_hinge(i,2) == 1
         res_DoFs = [res_DoFs, (2*i-1)];
     end
     if input_hinge(i,3) == 1
         res_DoFs = [res_DoFs, (2*i)];
     end
end
k_global = zeros(nnodes*2,nnodes*2,unknown);

for i = 1:nel
    type = input_pos(i,8);
    k_global(:,:,type) = k_global(:,:,type) + L(:,:,i)'*k_temp(:,:,i)*L(:,:,i);
end
u_observed = input_disp_b(:,2);

k_global(res_DoFs,:,:) =[];
k_global(:,res_DoFs,:) =[];

u_observed(res_DoFs) = [];
% Force vector
F = zeros(2*nnodes,1);
for i = 1:nnodes
    F(2*i-1) = input_load(i,2);
    F(2*i) = input_load(i,3);
end
F(res_DoFs) = [];
A = zeros((nnodes*2 - length(res_DoFs)),unknown);
for i = 1:unknown
    A(:,i) = k_global(:,:,i)*u_observed;
end

    
E_vec = zeros(4,1);
E_vec = (A'*A)\(A'*F);

E_out = cell(unknown+1,2);
E_out(1,:) = {'Element_Type','Youngs_Modulus(in Pa)'};
for i = 1:unknown
    E_out(i+1,2) = num2cell(ax_force(i));
end
E_out(1+1,1) = num2cell("Bottom");
E_out(2+1,1) = num2cell("Top");
E_out(3+1,1) = num2cell("Vertical");
E_out(4+1,1) = num2cell("Diagonal");
xlswrite('output_E.xlsx', E_out);

% E_vec = 

% Imposing Boundary condition


% % Different values of E
% E = zeros(1,4);
% coeff_E = zeros(nel,4);
% for cof = 1:length(E)
%     E = zeros(1,4);
%     E(cof) = 1;
%     for e = 1: nel
%         if input_pos(e,8) == 1
%             k(:,:,e) = (A*E(1)/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
%         end
%         if input_pos(e,8) == 2
%             k(:,:,e) = (A*E(2)/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
%         end
%         if input_pos(e,8) == 3
%             k(:,:,e) = (A*E(3)/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
%         end
%         if input_pos(e,8) == 4
%             k(:,:,e) = (A*E(4)/le(e)).*[1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0];
%         end
%     end
%     %% Assembly
%     K = zeros(2*nnodes,2*nnodes);
%     for e = 1:nel
%         K = K + ((L(:,:,e))'*(T(:,:,e))'*k(:,:,e)*T(:,:,e)*L(:,:,e));
%     end
% 
%     % Force vector
%     F = zeros(2*nnodes,1);
%     for i = 1:nnodes
%         F(2*i-1) = input_load(i,2);
%         F(2*i) = input_load(i,3);
%     end
% 
%     % Penalty Approach
%     K_hat = K;
%     F_hat = F;
%     k_mean = 0;
%     for i = 1:20
%         k_mean = k_mean + K(i,i); 
%     end
%     beta = k_mean * (10^7);
%     %hinge at position 1,2,12
%     for i = 1:nnodes
%         if input_hinge(i,2)== 1
%             K_hat((2*i-1),(2*i-1)) = beta + K(2*i-1);
%             F_hat(2*i-1) = beta*0;
%     %         K_hat((2*i-1),(2*i-1)) = 1;
%     %         F_hat = F_hat - (0.*K((2*i-1),:))';
%         end
%         if input_hinge(i,3)== 1
%             K_hat ((2*i),(2*i)) = beta + K(2*i);
%             F_hat(2*i) = beta*0;
%     %         K_hat((2*i),(2*i)) = 1;
%     %         F_hat =  F_hat - (0.*K((2*i),:))';
%         end
%     end
%     % Solution  
%     global_dof = K_hat\F_hat;
%     coeff_E(:,i) = global_dof; 
% end 
% 
% % E = coeff_E(1:4,:)\input_disp_b(1:4,2);
% 
% 
% % for i = 1:nnodes
% %     if input_hinge(i,2)== 1
% %         global_dof((2*i-1)) = 0;
% %     end
% %     if input_hinge(i,3)== 1
% %         global_dof((2*i)) = 0;
% %     end
% % end
% % %% Post Processing
% % local_dof = zeros(4,nel);
% % ax_force = zeros(nel,1);
% % for e = 1:nel
% %     local_dof(:,e) = T(:,:,e)*L(:,:,e)*global_dof;
% %     ax_force (e,1)= ((local_dof(3,e)-local_dof(1,e))*E*A)/(le(e));
% % end
%     
%     
% 
% 
% % Disp = cell(nnodes+1,3);
% % Disp(1,:) = {'Node_number','Displacement_x (in m)','Displacement_y (in m)'};
% % for i = 1:nnodes
% %     Disp(i+1,1) = num2cell(i);
% %     Disp(i+1,2) = num2cell(global_dof((2*i)-1));
% %     Disp(i+1,3) = num2cell(global_dof((2*i)));
% % end
% % Ax_force = cell(nel+1,2);
% % Ax_force(1,:) = {'Element_Number','Axial_Force (in N)'};
% % for i = 1:nel
% %     Ax_force(i+1,1) = num2cell(i);
% %     Ax_force(i+1,2) = num2cell(ax_force(i));
% % end
% % 
% % xlswrite('output_disp.xlsx', Disp);
% % xlswrite("output_ax_force", Ax_force);
% % %% Plotting
% % for e = 1:nel
% %     plot([input_pos(e,4),input_pos(e,6)],[input_pos(e,5),input_pos(e,7)],'k');
% %     axis equal;
% %     hold on;
% % end
% % for e = 1:nel
% %     i = input_pos(e,2);
% %     j = input_pos(e,3);
% %     plot([(input_pos(e,4)+global_dof(2*i-1)),(input_pos(e,6)+global_dof(2*j-1))],[(input_pos(e,5)+global_dof(2*i)),(input_pos(e,7)+global_dof(2*j))],'r');
% %     hold on;
% % end
% 
% 
% 
