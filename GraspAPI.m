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
                    obj.grasp_templates_folder = fullfile(cd,grasp_templates_folder);
                    obj.grasp_project_folder  = fullfile(cd,grasp_projects_folder,grasp_project_name);
                    obj.grasp_project_name = grasp_project_name;
            end
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
            gxp_template = obj.read_template('gxp.template');
            gxp_template = sprintf(gxp_template,strcat(obj.grasp_project_name,'.tor'),strcat(obj.grasp_project_name,'.tci'));
            output_file = fullfile(obj.grasp_project_folder,strcat(obj.grasp_project_name,'.gxp'));
            obj.create_write_file(output_file,gxp_template);
        end    
        function [] = create_global_coordinate_system(obj)
            %Create_frequency_axis Create a coordinate system with x,y,z
            %points and vectors
            %   Detailed explanation goes here
            gc_template_text = obj.read_template('coor_sys_global.template');
            obj.append_to_file(obj.grasp_tor_handler,gc_template_text);
        end
        function [] = create_coordinate_system(obj,name,origin,x_axis,y_axis,ref_coor)
            %Create_frequency_axis Create a coordinate system with x,y,z
            %points and vectors
            %   Detailed explanation goes here
            cs_template_text = obj.read_template('coor_sys.template');
            origin = obj.add_units(origin);
            cs_template_text = sprintf(cs_template_text,name,origin,x_axis,y_axis,ref_coor);
            obj.append_to_file(obj.grasp_tor_handler,cs_template_text);
        end
    end
    
    methods (Access = private)
        function [] = append_to_file(~,filehandler,text)
            fprintf(filehandler,"%s\n",text);
        end
        function [text] = read_template(obj,template_name)
            text = fileread(fullfile(obj.grasp_templates_folder,template_name));
        end
        function [filehandler] = check_create_file(obj,filename)
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
   end
end

