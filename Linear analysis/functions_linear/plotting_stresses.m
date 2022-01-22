function [ELEMENT]=plotting_stresses(ELEMENT,GEOMETRY,SF,type_SF)

xsi=[];
eta=[];
vect=[-1:0.5:1];
xsi(1,:)=vect;
eta(:,1)=flipud(vect');
for i=2:length(vect)
    xsi(i,:)=vect;
    eta(:,i)=flipud(vect');
end

for ii=1:GEOMETRY.N_elem
    if type_SF==4
        x=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).coord_nodes(1,1)+SF.N2(xsi,eta).*ELEMENT(ii).coord_nodes(2,1)+SF.N3(xsi,eta).*ELEMENT(ii).coord_nodes(3,1)+SF.N4(xsi,eta).*ELEMENT(ii).coord_nodes(4,1);
        y=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).coord_nodes(1,2)+SF.N2(xsi,eta).*ELEMENT(ii).coord_nodes(2,2)+SF.N3(xsi,eta).*ELEMENT(ii).coord_nodes(3,2)+SF.N4(xsi,eta).*ELEMENT(ii).coord_nodes(4,2);

        u=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).local_displ(1,1)+SF.N2(xsi,eta).*ELEMENT(ii).local_displ(2,1)+SF.N3(xsi,eta).*ELEMENT(ii).local_displ(3,1)+SF.N4(xsi,eta).*ELEMENT(ii).local_displ(4,1);
        v=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).local_displ(1,2)+SF.N2(xsi,eta).*ELEMENT(ii).local_displ(2,2)+SF.N3(xsi,eta).*ELEMENT(ii).local_displ(3,2)+SF.N4(xsi,eta).*ELEMENT(ii).local_displ(4,2);
    end
    if type_SF==8
        x=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).coord_nodes(1,1)+SF.N2(xsi,eta).*ELEMENT(ii).coord_nodes(2,1)+SF.N3(xsi,eta).*ELEMENT(ii).coord_nodes(3,1)+SF.N4(xsi,eta).*ELEMENT(ii).coord_nodes(4,1)+SF.N5(xsi,eta).*ELEMENT(ii).coord_nodes(5,1)+SF.N6(xsi,eta).*ELEMENT(ii).coord_nodes(6,1)+SF.N7(xsi,eta).*ELEMENT(ii).coord_nodes(7,1)+SF.N8(xsi,eta).*ELEMENT(ii).coord_nodes(8,1);
        y=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).coord_nodes(1,2)+SF.N2(xsi,eta).*ELEMENT(ii).coord_nodes(2,2)+SF.N3(xsi,eta).*ELEMENT(ii).coord_nodes(3,2)+SF.N4(xsi,eta).*ELEMENT(ii).coord_nodes(4,2)+SF.N5(xsi,eta).*ELEMENT(ii).coord_nodes(5,2)+SF.N6(xsi,eta).*ELEMENT(ii).coord_nodes(6,2)+SF.N7(xsi,eta).*ELEMENT(ii).coord_nodes(7,2)+SF.N8(xsi,eta).*ELEMENT(ii).coord_nodes(8,2);

        u=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).local_displ(1,1)+SF.N2(xsi,eta).*ELEMENT(ii).local_displ(2,1)+SF.N3(xsi,eta).*ELEMENT(ii).local_displ(3,1)+SF.N4(xsi,eta).*ELEMENT(ii).local_displ(4,1)+SF.N5(xsi,eta).*ELEMENT(ii).local_displ(5,1)+SF.N6(xsi,eta).*ELEMENT(ii).local_displ(6,1)+SF.N7(xsi,eta).*ELEMENT(ii).local_displ(7,1)+SF.N8(xsi,eta).*ELEMENT(ii).local_displ(8,1);
        v=@(xsi,eta) SF.N1(xsi,eta).*ELEMENT(ii).local_displ(1,2)+SF.N2(xsi,eta).*ELEMENT(ii).local_displ(2,2)+SF.N3(xsi,eta).*ELEMENT(ii).local_displ(3,2)+SF.N4(xsi,eta).*ELEMENT(ii).local_displ(4,2)+SF.N5(xsi,eta).*ELEMENT(ii).local_displ(5,2)+SF.N6(xsi,eta).*ELEMENT(ii).local_displ(6,2)+SF.N7(xsi,eta).*ELEMENT(ii).local_displ(7,2)+SF.N8(xsi,eta).*ELEMENT(ii).local_displ(8,2);
    end
    ELEMENT(i).x_disp=u;
    ELEMENT(i).y_disp=v;

%=============Plotting sigma_xy=======================================

    fprintf('Plotting sigma_xx for %d \n',ii)
    subplot(2,3,1)
    contourf(x(xsi,eta),y(xsi,eta),ELEMENT(ii).sigma_xx(xsi,eta)), hold on;
    c=colorbar;
    c.Label.String='MPa';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title(['\sigma_{xx}'])
    axis equal

%=============Plotting sigma_xy=======================================

    fprintf('Plotting sigma_yy for %d \n',ii)
    subplot(2,3,2)
    contourf(x(xsi,eta),y(xsi,eta),ELEMENT(ii).sigma_yy(xsi,eta)), hold on;
    c=colorbar;
    c.Label.String='MPa';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title(['\sigma_{yy}'])
    axis equal

%=============Plotting sigma_xy=======================================

    fprintf('Plotting sigma_xy for %d \n',ii)
    subplot(2,3,3)
    contourf(x(xsi,eta),y(xsi,eta),ELEMENT(ii).sigma_xy(xsi,eta)), hold on;
    c=colorbar;
    c.Label.String='MPa';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title(['\sigma_{xy}'])
    axis equal

%=============Plotting x displacements======================================= 
    
    fprintf('Plotting x displacements for %d \n',ii)
    subplot(2,3,4)
    contourf(x(xsi,eta),y(xsi,eta),u(xsi,eta)), hold on;
    c=colorbar;
    c.Label.String='mm';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title('u','fontsize',15,'interpreter','latex');
    axis equal
    
%=============Plotting y displacements======================================= 
    
    fprintf('Plotting x displacements for %d \n',ii)
    subplot(2,3,5)
    contourf(x(xsi,eta),y(xsi,eta),v(xsi,eta)), hold on;
    c=colorbar;
    c.Label.String='mm';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title('v','fontsize',15,'interpreter','latex');
    axis equal

%=============Plotting displacements magnitude======================================= 
    
    fprintf('Plotting displacements magnitude for %d \n',ii)
    subplot(2,3,6)
    contourf(x(xsi,eta),y(xsi,eta),sqrt(u(xsi,eta).^2+v(xsi,eta).^2)), hold on;
    c=colorbar;
    c.Label.String='mm';
    xlabel('x (mm)','fontsize',15,'interpreter','latex');
    ylabel('y (mm)','fontsize',15,'interpreter','latex');
    title('Disp','fontsize',15,'interpreter','latex');
    axis equal
    
end % END lement cycle 

end % END function 