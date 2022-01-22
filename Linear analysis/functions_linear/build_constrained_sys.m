function [K,F,Delete]=build_constrained_sys(GEOMETRY,K_glob_unc,F_glob_unc)

%=========Evaluate the constrained dof's================
Delete=[];    % vector that contains the constrained dof's
for i=1:length(GEOMETRY.spc(:,1))
    
    if GEOMETRY.spc(i,2)==1
        
        Delete(i)=(GEOMETRY.spc(i,1)-1)*2+1;
    end
    if GEOMETRY.spc(i,2)==2
       
        Delete(i)=(GEOMETRY.spc(i,1)-1)*2+2;
    end

end
%=======================================================

%==========Cancel out the rows and coloums related to constrained dof's====
K_glob_unc(Delete,:)=[];
K_glob_unc(:,Delete)=[];
F_glob_unc(Delete)=[];

K=K_glob_unc;
F=F_glob_unc;

end  % END function 
