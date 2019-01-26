function force = computeForceFieldSSD(Template, Source, displacement, gridObject)
x = gridObject.x;
y = gridObject.y;

Ux = displacement.x;
Uy = displacement.y;

dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

SSD = 0.0;

T = Template;
S = Source;

maximalBound = length(T);

% compute the unscaled SSD between template and source
for i = 2 : length(x)
    for j = 2: length(y)
        xpos = x(i); ypos = y(j);
        SSD = SSD + (1/2)*((T(xpos, ypos) - S(xpos, ypos)) * (T(xpos, ypos) - S(xpos, ypos)));
    end
end

% Sum Of Squared Difference
SSD = SSD / (length(x)-1) / (length(y) - 1);
disp(SSD);

for i = 2: length(x)
    for j = 2: length(y)
        xpos = x(i);
        ypos = y(j);
        
        % obtain the displaced points (x - u.x(x,ti)) and (y - u.y(y, ti))
        % which signifies x - u(x, ti) where x is treated as a vector of
        % the form x = [x1, x2]
        xpos_disp = ceil(x(i) - Ux(i, j));
        ypos_disp = ceil(y(j) - Uy(i, j));
        
        % add check making sure that U.x or U.y values dont push a point
        % outside boundaries 
        if(xpos_disp > 1 && ypos_disp > 1 && xpos_disp < maximalBound && ypos_disp < maximalBound)
            % coputer the intensity gradient
            dTdx = (T(xpos_disp+1, ypos_disp) - T(xpos_disp-1, ypos_disp)) / dx2;
            dTdy = (T(xpos_disp, ypos_disp+1) - T(xpos_disp, ypos_disp-1)) / dy2;

            % compute force field (-[T(x-u(x,t)) - S(x)]*Grad(T(x-U(x,t))))
            fx(i, j) = -(1/2)*(T(xpos_disp, ypos_disp) - S(xpos, ypos)) * dTdx;
            fy(i, j) = -(1/2)*(T(xpos_disp, ypos_disp) - S(xpos, ypos)) * dTdy;
        end
    end
end

force = struct();
force.x = fx;
force.y = fy;

end

