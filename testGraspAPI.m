clc;
clear;

%% Profile creation
baseTriang1 = 9; %cm
baseTriang2 = 7; %cm
htriang1 = 5; %cm
htriang2 = 8; %cm
broca = 0; %cm
ntriang = 8; 

triang_init = triangle_unit_broca(baseTriang1,baseTriang2,htriang1,htriang2,[0,0],[0,broca/2]);
triang_middle = triangle_unit_broca(baseTriang1,baseTriang2,htriang1,htriang2,[triang_init(end,1),0],[broca/2,broca/2]); 
for i = 1:(ntriang-3)
    triang_middle = [triang_middle;triangle_unit_broca(baseTriang1,baseTriang2,htriang1,htriang2,[triang_middle(end,1) 0],[broca/2,broca/2])];
end
triang_end = triangle_unit_broca(baseTriang1,baseTriang2,htriang1,htriang2,[triang_middle(end,1) 0],[broca/2,0]);
%triang = [0, 0; 9, 0];
def_profile = build_profile_from_form_xy([triang_init;triang_middle;triang_end],[1,1]);

plot(def_profile(:,1),def_profile(:,2));
%Adaptation to grasp format
grasp_profile = def_profile(1:(end-1),:)'./100;

%% Generate Grasp Diagram
[THin,PHin,e_th,e_ph] = leeFicheroASC("M1_32032.RFF.ASC");
escribeDiagramaGrasp(THin,PHin,e_th,e_ph,"sondamonofreq.pattern",10);
%% Grasp API
%Reflector parameters
reflector_diameter = 1;
focal_distance = 1.5;

%Init API object
grasp_handler = GraspAPI('GraspTemplates','GraspProjects','PruebaCoordenadas','m','GHz',true);
%Create coordinate systems
%Global coordinate system
grasp_handler.create_global_coordinate_system();
%Probe coordinate system
grasp_handler.create_coordinate_system("reflector_center_coor",[0,reflector_diameter/2,0],...
    [1,0,0],[0,1,0],"single_global_coor")
grasp_handler.create_coordinate_system_euler("probe_coor",[0,0,1.5],...
    [90 180-22.7 90],"single_global_coor")
%Create paraboloid and rim
grasp_handler.create_surface_paraboloid("reflector_paraboloid",focal_distance);
grasp_handler.create_tabulated_rim_xy("reflector_rim",grasp_profile,"sequential",[0,0]);
grasp_handler.create_reflector("single_reflector","reflector_center_coor",...
    "reflector_paraboloid","reflector_rim")
%Frequency List
%grasp_handler.create_frequency_range("frequencies_samples",10,120,12);
grasp_handler.create_frequency_list("frequencies_list",20);
%PO simulation
grasp_handler.create_po_analysis("single_po","frequencies_list","single_reflector","po_plus_ptd",...
    "single_global_coor");
%Output fields
grasp_handler.create_field_planar_cut("plane_cut_225","reflector_center_coor",2.25,[-0.5 0.5 201],[0 90 2],...
    "plane_cut_225.cut","frequencies_list");
grasp_handler.create_tabulated_feed("single_feed","frequencies_list","probe_coor","sonda32032.pattern",72,"far","off");
%Commands
grasp_handler.new_cmd_get_currents("single_po",["single_feed"],-80,"on",["plane_cut_225"]);
grasp_handler.new_cmd_get_field("plane_cut_225",["single_feed","single_po"]);

%Create GXP File for visualization
grasp_handler.create_gxp_file();
%Delete the object when is no longer used
grasp_handler.delete();


