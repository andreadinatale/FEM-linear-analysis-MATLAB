%% FE model 

clear all
close all
clc;

%% IMPORT GEOMETRY WITH MESH 

addpath("functions_linear")
addpath("Matlab FEM Input files")

%==========CHOOSE ANALYSIS YOU WANT TO PERFORM==========================
GEOMETRY=input_open_hole_biquad_44_els;     
Gauss_number=3;
Ampl_factor=200;
%=======================================================================
Size=size(GEOMETRY.elements);
type_SF=Size(2);
clear Size


x_matrix=GEOMETRY.nodes(:,1);
y_matrix=GEOMETRY.nodes(:,2);
z_matrix=zeros(size(GEOMETRY.nodes(:,1)));
plot(x_matrix,y_matrix,'k*'), grid on, hold on;
num_ele=size(GEOMETRY.elements);
title('Number of elements: ',num2str(num_ele(1)),'interpreter','latex');
xlabel('x (mm)','fontsize',15,'interpreter','latex');
ylabel('y (mm)','fontsize',15,'interpreter','latex');
axis equal

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
clear num_nodes num_ele
%% RESHAPING THE ELEMENTS AND THE NODES

[GEOMETRY]=reshaping(GEOMETRY,type_SF);

%% IMPORT SHAPE FUNCTIONS AND THEIR DERIVATIVES

[SF,dSF_xsi,dSF_eta]=SF(type_SF);

%% FEM MATRIX CREATION

[ELEMENT,K_glob_unc,F_glob_unc]=build_K_matrix(dSF_xsi,dSF_eta,GEOMETRY,Gauss_number,type_SF,MATERIAL);

%% CONSTRAINED SYSTEM

[K,F,Delete]=build_constrained_sys(GEOMETRY,K_glob_unc,F_glob_unc);

%% SOLVE THE LINEAR SYSTEM

[q_uncon_matrix]=FEM_solver(K,F,Delete,GEOMETRY);
clear K F K_glob_unc F_glob_unc

%% POST PROCESSING

fprintf('POST-PROCESSING \n\n')
fprintf('Amplificaticon factor used: %d \n',Ampl_factor)
plot(x_matrix+Ampl_factor*q_uncon_matrix(:,1),y_matrix+Ampl_factor*q_uncon_matrix(:,2),'bo','LineWidth',3)
legend('Nodes (Initial)','Nodes (Current)','Location','southwest')

Umax=max(abs(q_uncon_matrix(:,1)))
Vmax=max(abs(q_uncon_matrix(:,2)))

%% STRESS RECOVERY

fprintf('STRESS RECOVERY \n\n')
[ELEMENT]=stress_recovery(q_uncon_matrix,MATERIAL,ELEMENT,GEOMETRY,Gauss_number);

%% PLOTTING THE STRESS DISTRIBUTION 

[ELEMENT]=plotting_stresses(ELEMENT,GEOMETRY,SF,type_SF);




