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
        grasp_frequency_unit
    end
    
    methods
        function obj = GraspAPI(grasp_templates_folder,grasp_projects_folder,grasp_project_name,distance_unit,frequency_unit,overwrite)
            %GRASPAPI Construct an instance of this class
            %   Detailed explanation goes here
            switch(nargin)
                case {1,2}
                    error("Not enough input arguments");
                case 3
                    overwrite = false;
                    obj.grasp_distance_unit = 'm';
                    obj.grasp_frequency_unit = 'GHz';
                case 4
                    overwrite = false;
                    obj.grasp_distance_unit = distance_unit;
                    obj.grasp_frequency_unit = 'GHz';
                case 5
                    overwrite = false;
                    obj.grasp_distance_unit = distance_unit;
                    obj.grasp_frequency_unit = frequency_unit;
                otherwise
                    obj.grasp_distance_unit = distance_unit;
                    obj.grasp_frequency_unit = frequency_unit;
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
            %Link .tor to .tci
            obj.append_to_file(obj.grasp_tci_handler,sprintf("FILES READ ALL %s.tor",grasp_project_name));
        end
        function delete(obj)
            %Delete function is used for closing connection with files and
            %ssh
            obj.append_to_file(obj.grasp_tci_handler,"QUIT");
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
            origin = obj.add_distance_units(origin);
            cs_template_text = sprintf(cs_template_text,name,origin,x_axis,y_axis,ref_coor);
            obj.append_to_file(obj.grasp_tor_handler,cs_template_text);
        end
        function [] = create_coordinate_system_euler(obj,name,origin,euler_angles,ref_coor)
            %Create_frequency_axis Create a coordinate system with x,y,z
            %points and euler angle for the rotation of the coordinate
            %system
            cs_template_text = obj.read_template(fullfile("GeometricalObjects","CoordinateSystems",'coor_sys_euler.template'));
            origin = obj.add_distance_units(origin);
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
            vertex = obj.add_distance_units(vertex);
            focal_distance = obj.add_distance_units(focal_distance);
            surf_template_text = sprintf(surf_template_text,name,focal_distance,vertex);
            obj.append_to_file(obj.grasp_tor_handler,surf_template_text);
        end
        
        function [] = create_tabulated_rim_xy(obj,name,point_list,point_ordering,polar_origin,translation,rotation_angle,interpolation)
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
                    translation = [0,0];
                    rotation_angle = 0.0;
                    interpolation = "linear";                    
                case 6
                    rotation_angle = 0.0;
                    interpolation = "linear";
                case 7
                    interpolation = "linear";
            end
            rim_filename  = obj.create_rim_file(name,"%2.3f , %2.3f corner_point\n",point_list);
            polar_origin = obj.add_distance_units(polar_origin);
            translation = obj.add_distance_units(translation);
            tab_rim_template_text = obj.read_template(fullfile("GeometricalObjects","Rims","Tabulated Rims","tabulated_rim_xy.template"));
            tab_rim_template_text = sprintf(tab_rim_template_text,name,rim_filename,obj.grasp_distance_unit,...
                point_ordering,polar_origin,translation,rotation_angle,interpolation);
            obj.append_to_file(obj.grasp_tor_handler,tab_rim_template_text);
        end
        function [] = create_reflector(obj,name,coor_sys,surfaces,rim)
            %create_reflector create a reflector geometrical object
            surfaces_list = strjoin(strrep('ref(%s)','%s',surfaces),",");
            reflector_template_text = obj.read_template(fullfile("GeometricalObjects","Scatterer","reflector.template"));
            reflector_template_text = sprintf(reflector_template_text,name,coor_sys,surfaces_list,rim);
            obj.append_to_file(obj.grasp_tor_handler,reflector_template_text);
        end
        function [] = create_frequency_range(obj,name,min_freq,max_freq,n_freqs)
            freq_template_text = obj.read_template(fullfile("ElectricalObjects","Frequency","frequency_range.template"));
            freq_template_text = sprintf(freq_template_text,name,obj.add_frequency_units(min_freq),obj.add_frequency_units(max_freq),n_freqs);
            obj.append_to_file(obj.grasp_tor_handler,freq_template_text);
        end
        function [] = create_frequency_list(obj,name,list_freqs)
            freq_template_text = obj.read_template(fullfile("ElectricalObjects","Frequency","frequency_list.template"));
            list_freqs = strjoin(string(obj.add_frequency_units(list_freqs)),',');
            freq_template_text = sprintf(freq_template_text,name,list_freqs);
            obj.append_to_file(obj.grasp_tor_handler,freq_template_text);
        end
        function [] = create_po_analysis(obj,name,freq_ref,scatterer_ref,method,coor_sys,po_points,ptd_points)
            switch(nargin)
                case {1,2,3,4,5}
                    error("Not enough arguments")
                case 6
                    po_points = [0,0];
                    ptd_points = [-1,0];
                case 7
                    ptd_points = [-1,0];
            end
            po_analysis_text = obj.read_template(fullfile("ElectricalObjects","POAnalysis","po_single_face_scatterer.template"));
            [ncols, ~] = size(ptd_points);
            ptd_string = strings(1,ncols);
            for i = 1:ncols
                ptd_string(i) = sprintf("struct(edge:%d,ptd:%d)",ptd_points(i,:));
            end
            ptd_string = strjoin(ptd_string,",");
            po_analysis_text = sprintf(po_analysis_text,name,freq_ref,scatterer_ref,method,po_points,ptd_string,coor_sys,"po_ptd_currents.cur");
            obj.append_to_file(obj.grasp_tor_handler,po_analysis_text);
        end
        
        function [] = create_field_planar_cut(obj,name,coor_sys,near_dist,rho_range,phi_range,filename,freq_ref)
            cut_template_text = obj.read_template(fullfile("ElectricalObjects","FieldStorage","Cut","planar_cut.template"));
            rho_np = ceil(rho_range(3));
            phi_np = ceil(phi_range(3));
            near_dist = obj.add_distance_units(near_dist);
            cut_template_text  = sprintf(cut_template_text,name,coor_sys,near_dist,rho_range(1:2),...
                rho_np,obj.grasp_distance_unit,phi_range(1:2),phi_np,filename,freq_ref);
            obj.append_to_file(obj.grasp_tor_handler,cut_template_text);
        end
        
        function [] = create_field_planar_grid(obj,name,coor_sys,near_dist,x_range,y_range,polarization,filename,freq_ref)
            cut_template_text = obj.read_template(fullfile("ElectricalObjects","FieldStorage","Grid","planar_grid.template"));
            x_np = ceil(x_range(3));
            y_np = ceil(y_range(3));
            near_dist = obj.add_distance_units(near_dist);
            cut_template_text  = sprintf(cut_template_text,name,coor_sys,near_dist,x_range(1:2),...
                x_np,obj.grasp_distance_unit,y_range(1:2),y_np,polarization,filename,freq_ref);
            obj.append_to_file(obj.grasp_tor_handler,cut_template_text);
        end
        
        function [] = create_tabulated_feed(obj,name,freq_ref,coor_sys,filename_feed,n_cuts,near_far,far_forced)
            feed_template_text = obj.read_template(fullfile("ElectricalObjects","Feed","TabulatedFeed","tabulated_pattern.template"));
            feed_template_text = sprintf(feed_template_text,name,freq_ref,coor_sys,filename_feed,n_cuts,near_far,far_forced);
            copyfile(fullfile(cd,"GraspFeeds",filename_feed),fullfile(obj.grasp_project_folder,filename_feed))
            obj.append_to_file(obj.grasp_tor_handler,feed_template_text);
        end
        
        %Commands
        function [] = new_cmd_get_currents(obj,target,sources,field_accuracy,auto_convergence,convergence_references)
            cmd_template_text = obj.read_template(fullfile("Commands","get_currents.template"));
            sources = strjoin(strrep('ref(%s)','%s',sources),",");
            convergence_references = strjoin(strrep('ref(%s)','%s',convergence_references),",");
            cmd_template_text = sprintf(cmd_template_text,target,sources,field_accuracy,auto_convergence,convergence_references);
            obj.append_to_file(obj.grasp_tci_handler,cmd_template_text);
        end
       function [] = new_cmd_get_field(obj,target,sources)
            cmd_template_text = obj.read_template(fullfile("Commands","get_field.template"));
            sources = strjoin(strrep('ref(%s)','%s',sources),",");
            cmd_template_text = sprintf(cmd_template_text,target,sources);
            obj.append_to_file(obj.grasp_tci_handler,cmd_template_text);
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
        function [num_unit] = add_distance_units(obj,values)
            num_unit = obj.add_units(values,obj.grasp_distance_unit);
        end
        function [num_unit] = add_frequency_units(obj,values)
            num_unit = obj.add_units(values,obj.grasp_frequency_unit);
        end
        function [num_unit] = add_units(~,values,unit)
            num_unit = strings(size(values));
            for i = 1:length(values)
                num_unit(i) = sprintf("%f %s",values(i),unit);
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
        function [n_lines] = get_n_lines(filename)
            [status, cmdout] = system(sprintf('find /c /v "" %s',filename));
            if(status~=1)
                scanCell = textscan(cmdout,'%s %s %u');
                n_lines = scanCell{3};
            else
                n_lines = -1;
            end
        end
    end
    methods(Static)
        function [data_fields, x_axis] = read_cut(filename)
            %Get number of lines with windows command
            nlines = get_n_lines(filename);
            %Find First Header
            desc = readmatrix("plane_cut_200.txt",'Range',"A2:G2");
            initial_rad = desc(1);
            increment_rad = desc(2);
            n_rad = desc(3);
            %phi = desc(4);
            cmps = desc(7);
            phi_cuts = nlines/(n_rad+2);
            %create matrix and read thru phicuts
            data_matrix = zeros(phi_cuts,n_rad,cmps*2);
            for i = 1:phi_cuts
                range = sprintf("A%d:F%d",2*i+n_rad*(i-1)+1,(2+n_rad)*i);
                data_matrix(i,:,:) = readmatrix("plane_cut_200.txt",'Range',range);
            end
            x_axis = (0:(n_rad-1))*increment_rad + initial_rad;
            data_fields = complex(data_matrix(:,:,1:2:(end-1)),data_matrix(:,:,2:2:(end)));
        end
        function [data_fields, x_axis,y_axis] = read_grid(filename)
            datafile = fopen(filename,"r");
            %Find Header, dont care about bla bla
            header_offset = 1;
            while true
                tline = fgetl(datafile);
                header_offset = header_offset +1;
                if tline == "++++"
                    break
                end
            end
            %Check file type
            file_type = fgetl(datafile);
            if strtrim(file_type) ~= "1"
                error("This is not a GRASP grid file");
            end
            fclose(datafile);
            %Read component options
            range = sprintf("A%d:D%d",header_offset+1,header_offset+1);
            components = readmatrix(filename,'Range',range);
            n_components = components(3);
            %We do not care about the center
            %Find axes limits
            range = sprintf("A%d:D%d",header_offset+3,header_offset+3);
            axes = readmatrix(filename,'Range',range);
            %Find number of rows, columns
            range = sprintf("A%d:C%d",header_offset+4,header_offset+4);
            points = readmatrix(filename,'Range',range);
            %Create Axes
            x_axis = linspace(axes(1),axes(3),points(1));
            y_axis = linspace(axes(2),axes(4),points(2));
            if n_components == 2
                range = sprintf("A%d:D%d",header_offset+5,header_offset+5+points(1)*points(2));
            elseif n_components ==3
                range = sprintf("A%d:F%d",header_offset+5,header_offset+5+points(1)*points(2));
            else
                error("Unknown number of components")
            end
            data_matrix = readmatrix(filename,'Range',range);
            data_matrix = reshape(data_matrix,[points(1) points(2) n_components*2]);
            %Convert to complex
            data_fields = complex(data_matrix(:,:,1:2:(end-1)),data_matrix(:,:,2:2:(end)));
        end
    end
end

