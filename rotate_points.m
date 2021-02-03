function [points_rotated] = rotate_points(points,angle)
%rotate_points rotate an array of points by an angle
    rot_mat = [cosd(angle) -sind(angle);sind(angle) cosd(angle)];
    points_rotated = points*rot_mat;
end