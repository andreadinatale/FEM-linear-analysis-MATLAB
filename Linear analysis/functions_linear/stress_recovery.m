function [ELEMENT]=stress_recovery(q_uncon_matrix,MATERIAL,ELEMENT,GEOMETRY,Gauss_number)

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


for i=1:GEOMETRY.N_elem
    ELEMENT(i).local_displ=q_uncon_matrix(GEOMETRY.elements(i,:)',:);
    local_displ_vector=[];
    for k=1:length(ELEMENT(i).local_displ)
        local_displ_vector=[local_displ_vector;(ELEMENT(i).local_displ(k,:))'];
    end

%==============Calculating stresses at Gauss points=======================

    fprintf('Calculating stresses at Gauss points for %d element... \n',i);
    stress_Gauss=struct();
    for row=1:Gauss_number
        for column=1:Gauss_number
            sigma=MATERIAL.Q*ELEMENT(i).B(row,column).B*local_displ_vector;
            stress_Gauss(row,column).stress_Gauss=sigma;
        end
    end
%     ELEMENT(i).stress_Gauss=stress_Gauss;

%==============Calculating stress distributions on the element========================== 

    fprintf('Calculating stress distributions on the %d element... \n',i);
    sigma_xx_vect=[];
    sigma_yy_vect=[];
    sigma_xy_vect=[];
    Matrix=[];
    for row=1:Gauss_number
        for column=1:Gauss_number
                 sigma_xx_vect=[sigma_xx_vect; stress_Gauss(row,column).stress_Gauss(1)];
                 sigma_yy_vect=[sigma_yy_vect; stress_Gauss(row,column).stress_Gauss(2)];
                 sigma_xy_vect=[sigma_xy_vect; stress_Gauss(row,column).stress_Gauss(3)];
                 if Gauss_number==2
                     Matrix=[Matrix;1, x(column), x(row), x(column)*x(row)];
                 end
                 if Gauss_number==3
                     Matrix=[Matrix;1, x(column), x(row), x(column)*x(row), x(column)^2, x(row)^2, x(column)^2*x(row), x(column)*x(row)^2, x(column)^2*x(row)^2];
                 end
        end
    end
    a_xx=Matrix\sigma_xx_vect;
    a_yy=Matrix\sigma_yy_vect;
    a_xy=Matrix\sigma_xy_vect;
    if Gauss_number==2
        ELEMENT(i).sigma_xx=@(xsi,eta) a_xx(1)+a_xx(2).*xsi+a_xx(3).*eta+a_xx(4).*xsi.*eta;
        ELEMENT(i).sigma_yy=@(xsi,eta) a_yy(1)+a_yy(2).*xsi+a_yy(3).*eta+a_yy(4).*xsi.*eta;
        ELEMENT(i).sigma_xy=@(xsi,eta) a_xy(1)+a_xy(2).*xsi+a_xy(3).*eta+a_xy(4).*xsi.*eta;
    end
    if Gauss_number==3
        ELEMENT(i).sigma_xx=@(xsi,eta) a_xx(1)+a_xx(2).*xsi+a_xx(3).*eta+a_xx(4).*xsi.*eta+a_xx(5).*xsi.^2+a_xx(6).*eta.^2+a_xx(7).*xsi.^2.*eta+a_xx(8).*xsi.*eta.^2+a_xx(9).*xsi.^2.*eta.^2;
        ELEMENT(i).sigma_yy=@(xsi,eta) a_yy(1)+a_yy(2).*xsi+a_yy(3).*eta+a_yy(4).*xsi.*eta+a_yy(5).*xsi.^2+a_yy(6).*eta.^2+a_yy(7).*xsi.^2.*eta+a_yy(8).*xsi.*eta.^2+a_yy(9).*xsi.^2.*eta.^2;
        ELEMENT(i).sigma_xy=@(xsi,eta) a_xy(1)+a_xy(2).*xsi+a_xy(3).*eta+a_xy(4).*xsi.*eta+a_xy(5).*xsi.^2+a_xy(6).*eta.^2+a_xy(7).*xsi.^2.*eta+a_xy(8).*xsi.*eta.^2+a_xy(9).*xsi.^2.*eta.^2;
    end
  
end % END element cycle



end % END function 