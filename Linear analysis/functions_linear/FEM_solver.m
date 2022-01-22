function [q_uncon_matrix]=FEM_solver(K,F,Delete,GEOMETRY)

q=K\F;
q_uncon=zeros(2*GEOMETRY.N_nodes,1);
non_zero_pos=[1:2*GEOMETRY.N_nodes];
non_zero_pos(Delete)=[];
q_uncon(non_zero_pos)=q;
for i=1:GEOMETRY.N_nodes
     q_uncon_matrix(i,:)=[q_uncon(2*i-1) q_uncon(2*i)];
end

end