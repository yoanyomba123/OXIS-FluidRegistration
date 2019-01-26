function [perturbation, delta] = computePertubation(gridObject, U, V, tolerance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

perturbation = struct();
dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

Vx = V.x; Ux = U.x;
Vy = V.y; Uy = U.y;

% pad Arrays
Vx =padarray(Vx, [1,1]);  Vy = padarray(Vy, [1,1]);
Ux = padarray(Ux, [1,1]); Uy = padarray(Uy, [1,1]);

maxPerturb = 0.0;
perturbL2Norm = 0.0;
maxPertubation = 0.0;

for i = 2:length(Vx)-2
    for j= 2:length(Vy)-2
        dUx_dx = (Ux(i+1, j) - Ux(i-1, j))/dx2;
        dUx_dy = (Ux(i, j+1) - Ux(i, j-1))/dy2;

        
        dUy_dy = (Uy(i, j+1) - Uy(i, j-1))/dy2;
        dUy_dx = (Uy(i+1, j) - Uy(i-1, j))/dx2;

        perturbation.x(i, j) = Vx(i, j) - (Vx(i,j).* dUx_dx)  - (Vy(i, j) .* dUy_dx); %dUy_dx
        perturbation.y(i, j) = Vy(i, j) - (Vy(i, j).* dUy_dy) - (Vx(i,j) .* dUx_dy); %dUx_dy
        
        perturbL2Norm = sqrt(perturbation.x(i, j) .^2 +  perturbation.y(i, j) .^ 2);

        if perturbL2Norm > maxPertubation && perturbL2Norm ~= 0
            maxPertubation = perturbL2Norm;
        end
        
        if maxPertubation == 0.0
           maxPertubation = perturbL2Norm;
        end
    end
end

delta = tolerance.distanceTolerance / maxPertubation;
end

function mag = pertubationAt(x, y, U, V, Pertubation)
    % obtain pertubation at given point in image
    vi_x = V.x(x, y);
    vi_y = V.y(x, y);
    vi = [vi_x, vi_y];
    
    uxx = (U.x(x+1,y) - U.x(x-1, y));
    uyx = (U.y(x + 1) - U.y(x - 1, y));
    
    ux = [uxx, uyx] .* vi_x;
    
    uxy = (U.x(x, y+1) - U.x(x, y-1));
    uyy = (U.y(x, y+1) - U.y(x, y-1));
    
    uy = [uxy, uyy] .* vi_y;
    
    ue = (ux + uy) .* (1/2);
    
    ri = vi - ue;
    
    Pertubation.x(x,y) = ri(1);
    Pertubation.y(x,y) = ri(2);
    
    mag = norm(ri);
end