function [points_moved] = move_points(points,movement_vector)
%move_points move an array of points by a movement vector(x,y)
    movement_mat = ones(size(points)).*movement_vector;
    points_moved = points+movement_mat;
end