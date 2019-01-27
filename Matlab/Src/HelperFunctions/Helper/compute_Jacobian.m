function minimum = compute_Jacobian(gridObject, U, V, tolerance)
maxDeformation = tolerance.jacobianTolerance;

dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

Ux = U.x;
Uy = U.y;

% pad Arrays
Ux = padarray(Ux, [1,1]); Uy = padarray(Uy, [1,1]);
JacobianArray = struct();
minimum = 0.0;
count = 1;
for i = 2:length(Ux)-2
    for j= 2:length(Uy)-2
        dUx_dx = (Ux(i+1, j) - Ux(i-1, j))/dx2;
        dUx_dy = (Ux(i, j+1) - Ux(i, j-1))/dy2;
        dUy_dy = (Uy(i, j+1) - Uy(i, j-1))/dy2;
        dUy_dx = (Uy(i+1, j) - Uy(i-1, j))/dx2;
        
        jacobian = [dUx_dx, dUx_dy; dUy_dx, dUy_dy];
        jacobianDeterminant = det(jacobian);
        
        if count == 1
           minimum = jacobianDeterminant;
        end
        
        if jacobianDeterminant < minimum
           minimum = jacobianDeterminant;
        end
        
        count = count + 1;
    end
end


end

