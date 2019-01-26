function force = computeForceFieldJacMaps(Template, Source, displacement, gridObject)
% compute the force field by also taking into account jacobian maps in
% order to perserve topology once performing deformations;

x = gridObject.x;
y = gridObject.y;

U = displacement;
Ux = displacement.x;
Uy = displacement.y;

dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

SSD = 0.0;

T = Template;
S = Source;

lagrangeMult = 1;

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

for i = 2: length(U.x)
    for j = 2: length(U.y)
        xpos = x(i);
        ypos = y(j);
        
        % obtain the displaced points (x - u.x(x,ti)) and (y - u.y(y, ti))
        % which signifies x - u(x, ti) where x is treated as a vector of
        % the form x = [x1, x2]
        xpos_disp = ceil(real(x(i) - Ux(i, j)));
        ypos_disp = ceil(real(y(j) - Uy(i, j)));
        
        % add check making sure that U.x or U.y values dont push a point
        % outside boundaries 
        if(xpos_disp > 1 && ypos_disp > 1 && xpos_disp < maximalBound && ypos_disp < maximalBound)
            % coputer the intensity gradient
            dTdx = (T(xpos_disp+1, ypos_disp) - T(xpos_disp-1, ypos_disp)) / dx2;
            dTdy = (T(xpos_disp, ypos_disp+1) - T(xpos_disp, ypos_disp-1)) / dy2;
                
                if(i < length(U.x) && j < length(U.y))
                    dUxdx = (U.x(i+1, j) - U.x(i-1,j)) / dx2;
                    dUxdy = (U.x(i, j+1) - U.x(i,j-1)) / dy2;

                    dUydx = (U.y(i+1, j) - U.y(i-1,j)) / dx2;
                    dUydy = (U.y(i, j+1) - U.y(i,j-1)) / dy2;

                    jacDet = (dUxdx * dUydy) - (dUydx * dUydy);
                    if(jacDet == 0)
                        L = 0;
                    else
                        L = 1 + log(jacDet) - (1/jacDet);
                    end
                    % obtain discretized 2nd order derivative
                    dUxdxx = (U.x(i+1, j) - 2 * U.x(i,j) + U.x(i-1,j)) / (dx2^2);
                    dUxdyy = (U.x(i, j+1) - 2 * U.x(i,j) + U.x(i,j -1)) / (dy2^2);

                    dUydxx = (U.y(i+1, j) - 2 * U.y(i,j) + U.y(i-1,j)) / (dx2^2);
                    dUxdyy = (U.y(i, j+1) - 2 * U.y(i,j) + U.y(i,j-1)) / (dy2^2);

                    % partial of Uy W.R.T partial of x and partial of y
                    dUydxdy = (((U.y(i+1, j+1) - U.y(i+1, j-1)) .* L) - ((U.y(i-1, j+1) + U.y(i-1, j-1)) * L))  / (4 * dx * dy);
                    dUydydx = (((U.y(i+1, j+1) - U.y(i-1, j+1)) * L) - ((U.y(i+1, j-1) + U.y(i-1, j-1)) * L))  / (4 * dx * dy);


                    dUxdxdy = (((U.x(i+1, j+1) - U.x(i+1, j-1)) * L) - ((U.x(i-1, j+1) + U.x(i-1, j-1)) * L))  / (4 * dx * dy);
                    dUxdydx = (((U.x(i+1, j+1) - U.x(i-1, j+1)) * L) - ((U.x(i+1, j-1) + U.x(i-1, j-1)) * L))  / (4 * dx * dy);

                    jacMapx = -dUydxdy + dUydydx;
                    jacMapy = dUxdxdy - dUxdydx;

                    % compute force field (-[T(x-u(x,t)) - S(x)]*Grad(T(x-U(x,t))))
                    fx(i, j) = -(1/2)*(T(xpos_disp, ypos_disp) - S(xpos, ypos)) * dTdx - (lagrangeMult * jacMapx);
                    fy(i, j) = -(1/2)*(T(xpos_disp, ypos_disp) - S(xpos, ypos)) * dTdy - (lagrangeMult * jacMapx);
                end
            
        end
    end
end

force = struct();

fx = padarray(fx, [1,1], 0, "pre");

fy = padarray(fy, [1,1], 0, "pre");

force.x = fx;
force.y = fy;

end

