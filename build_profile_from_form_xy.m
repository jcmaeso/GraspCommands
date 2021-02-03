function [def_profile] = build_profile_from_form_xy(form,repetition_factor)
%build_profile_from_form This function builds a profile from an input form
%with a repetition factor
%   Detailed explanation goes here
    profile_top = build_profile_side(form,repetition_factor(1));
    movement_dist_x = profile_top(end,1)-profile_top(1,1);
    
    profile_right = move_points(rotate_points(build_profile_side(form,repetition_factor(2)),90),[movement_dist_x,0]);
    movement_dist_y = profile_right(end,2)-profile_right(1,2);
    profile_bottom = move_points(rotate_points(profile_top,180),[movement_dist_x,movement_dist_y]);
    profile_left = move_points(rotate_points(profile_right,180),[movement_dist_x,movement_dist_y]);
    %We need to omit common points
    def_profile = move_points([profile_top;profile_right;profile_bottom;profile_left],[-movement_dist_x/2,-movement_dist_y/2]);
end


function [profile] = build_profile_side(form,repetitions)
%move_triangle Move Triangle by XY vector
    offset = ones(size(form)).*form(end,:);
    profile = form;
    for i = 1:repetitions-1
        new_form = form+i*offset;
        profile = [profile;new_form(2:end,:)];
    end
end


