clc;
clear;
%Init API object
grasp_handler = GraspAPI('GraspTemplates','GraspProjects','PruebaCoordenadas','m',true);
%Create coordinate systems
%Global coordinate system
grasp_handler.create_global_coordinate_system();
%Probe coordinate system
grasp_handler.create_coordinate_system("reflector_centre_coor",[0,0,0],...
    [0.92307692407069, 0,0.384615382230345],[0,0,0],"single_global_coor")
grasp_handler.create_coordinate_system_euler("probe_coor",[0,0,1.5],...
    [0 180-22.7 180],"single_global_coor")
%Create GXP File for visualization
grasp_handler.create_gxp_file();
%Delete the object when is no longer used
grasp_handler.delete();

