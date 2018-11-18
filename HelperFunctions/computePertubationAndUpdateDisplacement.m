function U = computePertubationAndUpdateDisplacement(gridObject, U, V, tolerance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

maxDeformation = tolerance.deformationTolerance;

perturbation = struct();
dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

Vx = V.x; Ux = U.x;
Vy = V.y; Uy = U.y;

% pad Arrays
Vx = padarray(Vx, [1,1]); Vy = padarray(Vy, [1,1]);
Ux = padarray(Ux, [1,1]); Uy = padarray(Uy, [1,1]);

maxPerturb = 0.0;

for i = 2:length(Vx)-2
    for j= 2:length(Vy)-2
        dUx_dx = (Ux(i+1, j) - Ux(i-1, j))/dx2;
        
        dUy_dy = (Uy(i, j+1) - Uy(i, j-1))/dy2;

        perturbation.x(i, j) = Vx(i, j) - Vx(i,j).*dUx_dx - Vy(i, j).*dUy_dy; 
        perturbation.y(i, j) = Vy(i, j) - Vx(i,j).*dUx_dx - Vy(i, j).*dUy_dy; 
        
        perturbL2Norm = sqrt(perturbation.x(i, j) .^2 +  perturbation.y(i, j) .^ 2);
        
        if(perturbL2Norm <= maxDeformation && perturbL2Norm ~= 0)
            maxPerturb = perturbL2Norm;
        else
            maxPerturb = maxDeformation;
        end
    end
end

dt = maxPerturb;
U.x = U.x + perturbation.x .* dt;
U.y = U.y + perturbation.y .* dt;

end

