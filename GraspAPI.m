classdef GraspAPI
    %GRASPAPI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        grasp_templates_folder
        grasp_project_folder
        grasp_project_name
        grasp_tci_handler
        grasp_tor_handler
        grasp_distance_unit
    end
    
    methods
        function obj = GraspAPI(grasp_templates_folder,grasp_projects_folder,grasp_project_name,distance_unit,overwrite)
            %GRASPAPI Construct an instance of this class
            %   Detailed explanation goes here
            switch(nargin)
                case {1,2}
                    error("Not enough input arguments");
                case 3
                    overwrite = false;
                    obj.grasp_distance_unit = 'm';
                case 4
                    overwrite = false;
                    obj.grasp_distance_unit = distance_unit;
                otherwise
                    obj.grasp_distance_unit = distance_unit;
            end
                    obj.grasp_templates_folder = fullfile(cd,grasp_templates_folder);
                    obj.grasp_project_folder  = fullfile(cd,grasp_projects_folder,grasp_project_name);
                    obj.grasp_project_name = grasp_project_name;
            %Check if project folder exists, if not create it
            if ~exist(grasp_projects_folder, 'dir')
               mkdir(grasp_projects_folder);
            end
            if ~exist(obj.grasp_project_folder, 'dir')
                mkdir(obj.grasp_project_folder)
            end
            if overwrite
                obj.grasp_tci_handler = obj.overwrite_create_file(fullfile(obj.grasp_project_folder,strcat(grasp_project_name,'.tci')));
                obj.grasp_tor_handler = obj.overwrite_create_file(fullfile(obj.grasp_project_folder,strcat(grasp_project_name,'.tor')));
            else
                %Check if tci/tor file exists, if not create it
                obj.grasp_tci_handler = obj.check_create_file(fullfile(obj.grasp_project_folder,strcat(grasp_project_name,'.tci')));
                obj.grasp_tor_handler = obj.check_create_file(fullfile(obj.grasp_project_folder,strcat(grasp_project_name,'.tor')));
            end
        end
        function delete(obj)
            %Delete function is used for closing connection with files and
            %ssh
            fclose(obj.grasp_tci_handler);
            fclose(obj.grasp_tor_handler);
        end
        function [] = create_gxp_file(obj)
            gxp_template = obj.read_template(fullfile("Files","gxp.template"));
            gxp_template = sprintf(gxp_template,strcat(obj.grasp_project_name,'.tor'),strcat(obj.grasp_project_name,'.tci'));
            output_file = fullfile(obj.grasp_project_folder,strcat(obj.grasp_project_name,'.gxp'));
            obj.create_write_file(output_file,gxp_template);
        end
        
        function [] = create_global_coordinate_system(obj)
            %Create_frequency_axis Create a global coordinate system in
            %0,0,0
            %   Create basic coordinate system by initilizating an empty
            %   one used as reference for other primary coordinate systems
            gc_template_text = obj.read_template(fullfile("GeometricalObjects","CoordinateSystems",'coor_sys_global.template'));
            obj.append_to_file(obj.grasp_tor_handler,gc_template_text);
        end
        function [] = create_coordinate_system(obj,name,origin,x_axis,y_axis,ref_coor)
            %Create_frequency_axis Create a coordinate system with x,y,z
            %points and direction vectors for x and y axis
            cs_template_text = obj.read_template(fullfile("GeometricalObjects","CoordinateSystems",'coor_sys.template'));
            origin = obj.add_units(origin);
            cs_template_text = sprintf(cs_template_text,name,origin,x_axis,y_axis,ref_coor);
            obj.append_to_file(obj.grasp_tor_handler,cs_template_text);
        end
        function [] = create_coordinate_system_euler(obj,name,origin,euler_angles,ref_coor)
            %Create_frequency_axis Create a coordinate system with x,y,z
            %points and euler angle for the rotation of the coordinate
            %system
            cs_template_text = obj.read_template(fullfile("GeometricalObjects","CoordinateSystems",'coor_sys_euler.template'));
            origin = obj.add_units(origin);
            cs_template_text = sprintf(cs_template_text,name,origin,euler_angles,ref_coor);
            obj.append_to_file(obj.grasp_tor_handler,cs_template_text);
        end
        
        function [] = create_surface_paraboloid(obj,name,focal_distance,vertex)
            %create_surface_paraboloid creates a grasp paraboloid
            %geometrical object
            %----Input-----------
            %name: grasp name object (string)
            %focal_distance: paraboloid focal distance for graps object(decimal)
            %vertex: paraboloid vertex position, if not provided defaulted
            %to [0,0,0]
            switch(nargin)
                case {1,2}
                    error("Not enough input arguments");
                case 3
                    vertex = [0,0,0];
            end
            surf_template_text = obj.read_template(fullfile("GeometricalObjects","Surfaces",'paraboloid.template'));
            vertex = obj.add_units(vertex);
            focal_distance = obj.add_units(focal_distance);
            surf_template_text = sprintf(surf_template_text,name,focal_distance,vertex);
            obj.append_to_file(obj.grasp_tor_handler,surf_template_text);
        end
        
        function [] = create_tabulated_rim_xy(obj,name,point_list,point_ordering,polar_origin,rotation_angle,interpolation)
            %create_tabulated_rim_xy tabulated rim from a list of points
            %----Input-----------
            %name: grasp name object (string) and .rim file
            %point_list: list of points x-y format 2xN matrix
            %point_ordering: "sorted" or "sequential"
            %polar_origin: x-y center point
            %rotation_angle: rotation of points list (default 0,0)
            %interpolation: "linear"(default) or "spline"
            switch(nargin)
                case {1,2,3,4}
                    error("Not enough arguments")
                case 5
                    rotation_angle = 0.0;
                    interpolation = "linear";
                case 6
                    interpolation = "linear";
            end
            rim_filename  = obj.create_rim_file(name,"%2.3f , %2.3f corner_point\n",point_list);
            polar_origin = obj.add_units(polar_origin);
            tab_rim_template_text = obj.read_template(fullfile("GeometricalObjects","Rims","Tabulated Rims","tabulated_rim_xy.template"));
            tab_rim_template_text = sprintf(tab_rim_template_text,name,rim_filename,obj.grasp_distance_unit,...
                point_ordering,polar_origin,rotation_angle,interpolation);
            obj.append_to_file(obj.grasp_tor_handler,tab_rim_template_text);
        end
        function [] = create_reflector(obj,name,coor_sys,surfaces,rim)
            %create_reflector create a reflector geometrical object
            surfaces_list = strjoin(strrep('ref(%s)','%s',surfaces),",");
            reflector_template_text = obj.read_template(fullfile("GeometricalObjects","Scatterer","reflector.template"));
            reflector_template_text = sprintf(reflector_template_text,name,coor_sys,surfaces_list,rim);
            obj.append_to_file(obj.grasp_tor_handler,reflector_template_text);
        end
    end
    
    methods (Access = private)
        function [] = append_to_file(~,filehandler,text)
            fprintf(filehandler,"%s\n",text);
        end
        function [text] = read_template(obj,template_name)
            text = fileread(fullfile(obj.grasp_templates_folder,template_name));
        end
        function [filehandler] = check_create_file(~,filename)
            if exist(filename, 'file')
                filehandler = fopen(filename,'a');
            else
                filehandler = fopen(filename,'w');
            end
        end
        function [filehandler] = overwrite_create_file(~,filename)
            filehandler = fopen(filename,'w');
        end
        function [] = create_write_file(~,filename,text)
            filehandler = fopen(filename,'w');
            fprintf(filehandler,"%s\n",text);
            fclose(filehandler);
        end
        function [num_unit] = add_units(obj,values)
            num_unit = strings(size(values));
            for i = 1:length(values)
                num_unit(i) = sprintf("%f %s",values(i),obj.grasp_distance_unit);
            end
        end
        function [rim_filename] = create_rim_file(obj,name,line_template,list)
            %create_rim_file rim free format file from line template and
            %list of points generated by obj name and project name to
            %generate unicity
            rim_filename = strcat(obj.grasp_project_name,'_',name,'.rim');
            rim_file = obj.overwrite_create_file(fullfile(obj.grasp_project_folder,rim_filename));
            %Header
            fprintf(rim_file,"#X - Y\n");
            %List of points
            fprintf(rim_file,line_template,list);
            fclose(rim_file);
        end
   end
end

