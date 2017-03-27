function [x_coords, y_coords] = cylinderInterior(n,m,x0,y0,radius)
% Finds the x and y coordinates of all of the grid points in the interior
% of the cylinder and outputs two arrays: one for the x-coordinates (m) 
% and one for the y-coordinates (n).
x_coords = [];
y_coords = [];
for b = 1:n
    for a = 1:m
        dist = sqrt((a-x0)^2 + (b-y0)^2);
        if dist <= radius-1
            x_coords = cat(1,x_coords,a);
            y_coords = cat(1,y_coords,b);
        end
    end
end

end

