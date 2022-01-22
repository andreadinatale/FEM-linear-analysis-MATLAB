function [ELEMENT,K_glob_unc,F_glob_unc]=build_K_matrix(dSF_xsi,dSF_eta,GEOMETRY,Gauss_number,type_SF,MATERIAL)

%==============Cell array of shape functions and their derivatives=========
if type_SF==8
% vect_N={SF.N1 SF.N2 SF.N3 SF.N4 SF.N5 SF.N6 SF.N7 SF.N8};
vect_dN_xsi={dSF_xsi.dN1_xsi dSF_xsi.dN2_xsi dSF_xsi.dN3_xsi dSF_xsi.dN4_xsi dSF_xsi.dN5_xsi dSF_xsi.dN6_xsi dSF_xsi.dN7_xsi dSF_xsi.dN8_xsi};
vect_dN_eta={dSF_eta.dN1_eta dSF_eta.dN2_eta dSF_eta.dN3_eta dSF_eta.dN4_eta dSF_eta.dN5_eta dSF_eta.dN6_eta dSF_eta.dN7_eta dSF_eta.dN8_eta};
end
if type_SF==4
% vect_N={SF.N1 SF.N2 SF.N3 SF.N4};
vect_dN_xsi={dSF_xsi.dN1_xsi dSF_xsi.dN2_xsi dSF_xsi.dN3_xsi dSF_xsi.dN4_xsi};
vect_dN_eta={dSF_eta.dN1_eta dSF_eta.dN2_eta dSF_eta.dN3_eta dSF_eta.dN4_eta}; 
end
%==============Vector of Gauss point==============

if Gauss_number==3
x=GEOMETRY.int_rule.three_point.x;
end
if Gauss_number==2
x=GEOMETRY.int_rule.two_point.x;
end
if Gauss_number==1
x=GEOMETRY.int_rule.one_point.x;
end

%==============Vector of Gauss weigth==============

if Gauss_number==3
w=GEOMETRY.int_rule.three_point.w;
end
if Gauss_number==2
w=GEOMETRY.int_rule.two_point.w;
end
if Gauss_number==1
w=GEOMETRY.int_rule.one_point.w;
end
%======================================================================================================

K_glob_unc=zeros(2*GEOMETRY.N_nodes,2*GEOMETRY.N_nodes);
%=============Define ELEMENT struct===================

ELEMENT=struct();
for i=1:GEOMETRY.N_elem

     el=GEOMETRY.nodes(GEOMETRY.elements(i,:),:);
     ELEMENT(i).coord_nodes=el;

%=============Define Jacobian J at each Gauss point ===================
     
     fprintf('Calculating J at Gauss points for %d element... \n',i);
     coord_nodes=el';
     J_struct=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
              if type_SF==8
                  grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row)) vect_dN_xsi{5}(x(column),x(row)) vect_dN_xsi{6}(x(column),x(row)) vect_dN_xsi{7}(x(column),x(row)) vect_dN_xsi{8}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row)) vect_dN_eta{5}(x(column),x(row)) vect_dN_eta{6}(x(column),x(row)) vect_dN_eta{7}(x(column),x(row)) vect_dN_eta{8}(x(column),x(row))]'];
                  J_struct(row,column).J=coord_nodes*grad_xsi_eta;
              end
              if type_SF==4
                  grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row))]'];
                  J_struct(row,column).J=coord_nodes*grad_xsi_eta;
              end
         end
     end
     ELEMENT(i).J=J_struct;

%=============Define gradient in physical domain===================

     fprintf('Calculating gradient in physical domain at Gauss points for %d element... \n',i);
     GradN=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
             J=J_struct(row,column).J;
             if type_SF==8
                 grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row)) vect_dN_xsi{5}(x(column),x(row)) vect_dN_xsi{6}(x(column),x(row)) vect_dN_xsi{7}(x(column),x(row)) vect_dN_xsi{8}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row)) vect_dN_eta{5}(x(column),x(row)) vect_dN_eta{6}(x(column),x(row)) vect_dN_eta{7}(x(column),x(row)) vect_dN_eta{8}(x(column),x(row))]'];
                 GradN(row,column).GradN=J'\grad_xsi_eta';
             end
             if type_SF==4
                 grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row))]'];
                 GradN(row,column).GradN=J'\grad_xsi_eta';
             end
         end
     end
%      ELEMENT(i).GradN=GradN;

%=============Define F matrix=====================================     

%      fprintf('Calculating F matrix at Gauss points for %d element... \n',i);
%      F=struct();
%      for row=1:Gauss_number
%          for column=1:Gauss_number   
%              F(row,column).F=coord_nodes*(ELEMENT(i).GradN(row,column).GradN)';
%          end
%      end
%      ELEMENT(i).F=F;

%=============Define B matrix=====================================

     fprintf('Calculating B matrix at Gauss points for %d element... \n',i);
     B=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
             B_vect=[];
             for j=1:type_SF
                 grad=GradN(row,column).GradN(:,j);  
                 B_alpha=[grad(1) 0;0 grad(2);grad(2) grad(1)];
                 B_vect(:,[2*j-1 2*j])=B_alpha;
             end
             B(row,column).B=B_vect;
         end
     end
     ELEMENT(i).B=B;

%=============Define K_ele matrix=====================================     
             
     fprintf('Calculating K_ele matrix for %d element... \n',i);
     K_ele=zeros(2*type_SF,2*type_SF);
     for row=1:Gauss_number
         for column=1:Gauss_number
             B=ELEMENT(i).B(row,column).B;
             C=MATERIAL.Q;
             J=J_struct(row,column).J;
             K_ele=K_ele+(B'*C*B).*(det(J)*w(row)*w(column));
         end
     end
     ELEMENT(i).K_ele=K_ele;

%=============Define K matrix=====================================
     
     fprintf('Calculating K matrix for %d element... \n\n',i);
     K_unc=GEOMETRY.K_unc;
     K_unc(GEOMETRY.ptrs(i,:),GEOMETRY.ptrs(i,:))=K_ele;
     K_glob_unc=K_glob_unc+K_unc;

end  % END element cycle

%=============Define F vector=====================================

F_glob_unc=zeros(2*GEOMETRY.N_nodes,1);
for i=1:length(GEOMETRY.load(:,1))
F_glob_unc((GEOMETRY.load(i,1)-1)*2+1)=GEOMETRY.load(i,3);
end

end    % END function
