function [new_triangle,or_triangle] = triangle_unit_broca(triangle1_base,triangle2_base,triangle1_heigth,triangle2_heigth,origin,broca)
%move_triangle Move Triangle by XY vector
new_triangle = [0,0;...
    triangle1_base/2-triangle2_base/2,triangle1_heigth;...
    triangle1_base/2,triangle2_heigth+triangle1_heigth;...
    triangle1_base/2+triangle2_base/2,triangle1_heigth;
    triangle1_base,0];
or_triangle = new_triangle;

m_inicio = (new_triangle(2,2)-new_triangle(1,2))/(new_triangle(2,1)-new_triangle(1,1));
eq_recta_inicio = @(x) m_inicio*(x-new_triangle(2,1))+new_triangle(2,2);
cut_inicio = [broca(1) eq_recta_inicio(broca(1))];
beginning = [0 eq_recta_inicio(broca(1))];

m_fin = (new_triangle(end,2)-new_triangle(end-1,2))/(new_triangle(end,1)-new_triangle(end-1,1));
eq_recta_fin = @(x) m_fin*(x-new_triangle(end,1))+new_triangle(end,2);
cut_fin = [triangle1_base-broca(2) eq_recta_fin(triangle1_base-broca(2))];
ending = [triangle1_base eq_recta_fin(triangle1_base-broca(2))];
new_triangle = [beginning;cut_inicio;new_triangle(2:end-1,:);cut_fin;ending];
new_triangle(:,1)= new_triangle(:,1)+origin(1);
new_triangle(:,2)= new_triangle(:,2)+origin(2);
or_triangle(:,1) = or_triangle(:,1)+origin(1);
or_triangle(:,2) = or_triangle(:,2)+origin(2);
end