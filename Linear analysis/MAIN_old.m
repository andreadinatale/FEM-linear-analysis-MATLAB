%% FE model 

clear all
close all
clc;

%% IMPORT GEOMETRY WITH MESH 

addpath("my_functions")
addpath("Matlab FEM Input files")

%===============UPLAOD THE GEOMETRY JUST FOR PLOTTING THE ELEMENTS===============
INPUT_44 = input_open_hole_bilinear_44_els;    
INPUT_570 = input_open_hole_bilinear_570_els;
INPUT_2400 = input_open_hole_bilinear_2400_els;
%================================================================================

%===============UPLOAD ALL GEOMETRIES============================= 
GEOMETRY_S4_44 = input_open_hole_bilinear_44_els; 
GEOMETRY_S4_570 = input_open_hole_bilinear_570_els; 
GEOMETRY_S4_2400 = input_open_hole_bilinear_2400_els; 
GEOMETRY_S8_44 = input_open_hole_biquad_44_els; 
GEOMETRY_S8_570 = input_open_hole_biquad_570_els; 
GEOMETRY_S8_2400 = input_open_hole_biquad_2400; 
%=========================================================================

%==========CHOOSE NUMBER OF ELEMENTS (just for plotting)==========================
INPUT=INPUT_44;
%===========================================================================

x_matrix_S4=INPUT.nodes(:,1);
y_matrix_S4=INPUT.nodes(:,2);
z_matrix_S4=zeros(size(INPUT.nodes(:,1)));
figure('Name','Geometry','NumberTitle','off');
plot(x_matrix_S4,y_matrix_S4,'ko','LineWidth',5), grid on, hold on;
xlabel('x (mm)','fontsize',15,'interpreter','latex');
ylabel('y (mm)','fontsize',15,'interpreter','latex');
axis equal

%==========CHOOSE ANALYSIS YOU WANT TO PERFORM==========================
GEOMETRY=GEOMETRY_S4_44;      
type_SF=4;
Gauss_number=2;
Ampl_factor=100;
%=======================================================================

x_matrix=GEOMETRY.nodes(:,1);
y_matrix=GEOMETRY.nodes(:,2);
z_matrix=zeros(size(GEOMETRY.nodes(:,1)));
plot(x_matrix,y_matrix,'r*')
num_ele=size(GEOMETRY.elements);
title('Number of elements: ',num2str(num_ele(1)),'interpreter','latex');

%=============Define MATERIAL struct====================================
MATERIAL=struct();
MATERIAL.E=GEOMETRY.E;
MATERIAL.nu=GEOMETRY.nu;
MATERIAL.G=MATERIAL.E/(2*(1+MATERIAL.nu));
MATERIAL.t=GEOMETRY.t;
MATERIAL.Q=(MATERIAL.E*MATERIAL.t/(1-MATERIAL.nu^2)).*[1 MATERIAL.nu 0
    MATERIAL.nu 1 0
    0 0 (1-MATERIAL.nu)/2];
%=======================================================================

%===================CREATE GAUSS POINTS=================================
GEOMETRY.int_rule=struct();
GEOMETRY.int_rule.one_point=struct();
GEOMETRY.int_rule.one_point.x=[0];
GEOMETRY.int_rule.one_point.w=[2];
GEOMETRY.int_rule.two_point=struct();
GEOMETRY.int_rule.two_point.x=[-1/sqrt(3),1/sqrt(3)];
GEOMETRY.int_rule.two_point.w=[1,1];
GEOMETRY.int_rule.three_point=struct();
GEOMETRY.int_rule.three_point.x=[-sqrt(0.6),0,sqrt(0.6)];
GEOMETRY.int_rule.three_point.w=[5/9,8/9,5/9];
%======================================================================

GEOMETRY.N_elem=num_ele(1);
num_nodes=size(GEOMETRY.nodes);
GEOMETRY.N_nodes=num_nodes(1);

%% RESHAPING THE ELEMENTS AND THE NODES

[GEOMETRY]=reshaping(GEOMETRY,type_SF);

%% IMPORT SHAPE FUNCTIONS AND THEIR DERIVATIVES

[SF,dSF_xsi,dSF_eta]=SF(type_SF);

%% FEM MATRIX CREATION

[ELEMENT]=Fem_matrix_creation(SF,dSF_xsi,dSF_eta,GEOMETRY,Gauss_number,type_SF);

%% FEM MATRIX ASSEMBLY

[ELEMENT,K,F,Delete]=Fem_matrix_assembly(ELEMENT,MATERIAL,GEOMETRY,Gauss_number,type_SF);

%% FEM SOLVER

[q_uncon_matrix,q_uncon,q]=FEM_solver(K,F,Delete,GEOMETRY);

%% POST PROCESSING

fprintf('Amplificaticon factor used: %d \n',Ampl_factor)
plot(x_matrix+Ampl_factor*q_uncon_matrix(:,1),y_matrix+Ampl_factor*q_uncon_matrix(:,2),'b*')
legend('Element corners','Nodes (Initial)','Nodes (Current)')

Umax=max(abs(q_uncon_matrix(:,1)))
Vmax=max(abs(q_uncon_matrix(:,2)))

%% STRESS RECOVERY

for i=1:GEOMETRY.N_elem
    ELEMENT(i).local_displ=ELEMENT(i).Omega*q_uncon;
    stress_Gauss=struct();
    for row=1:Gauss_number
        for column=1:Gauss_number
            sigma=MATERIAL.Q*ELEMENT(i).B(row,column).B*ELEMENT(i).local_displ;
            stress_Gauss(row,column).stress_Gauss=sigma;
        end
    end
    ELEMENT(i).stress_Gauss=stress_Gauss;
    
end













