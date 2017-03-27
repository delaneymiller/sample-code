function [array] = reassignCylinder(array,value,x_coords,y_coords)
% Assigns specified value to all points in the cylinder in specified array.
% For example, if we want psi = 0 everywhere in the cylinder, we would
% say: psi = reassignCylinder(psi,0,x,y);
points = length(y_coords);
for a = 1:points
    xa = x_coords(a);
    ya = y_coords(a);
    array(ya,xa) = value;
end
end

