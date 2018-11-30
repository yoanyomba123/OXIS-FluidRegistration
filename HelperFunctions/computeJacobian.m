function regridBool = computeJacobian(gridObject, U, maxPertubation, perturbation, tolerance)
regridBool = "false";
dx = gridObject.dx;
dy = gridObject.dy;

dx2 = dx * 2.0;
dy2 = dy * 2.0;
% process the jacobian
jacobianObject = struct();
jacobianObject.referenceFieldEntityX = gridObject.grid.x - U.x - (maxPertubation .* perturbation.x);
jacobianObject.referenceFieldEntityY = gridObject.grid.y - U.y - (maxPertubation .* perturbation.y);

entityX = jacobianObject.referenceFieldEntityX;
entityY = jacobianObject.referenceFieldEntityY;

entityX = padarray(entityX, [1,1]); entityY = padarray(entityY, [1,1]);
jacobianAvg = 0.0;
for i = 2:length(entityX)-2
    for j= 2:length(entityX)-2     
        % process the jacobian
        dEntx_dx = (entityX(i+1, j) - entityX(i-1, j))/dx2;
        dEntx_dy = (entityX(i, j+1) - entityX(i, j-1))/dy2;

        dEnty_dy = (entityY(i, j+1) - entityY(i, j-1))/dy2;
        dEnty_dx = (entityY(i+1, j) - entityY(i-1, j))/dx2;
        
        Jacobian = [dEntx_dx, dEntx_dy; dEnty_dx, dEnty_dy];
        jacobianAvg = jacobianAvg + abs(det(Jacobian));
       
    end
end
jacobianAvg = (jacobianAvg/ length( 2:length(entityX)-2));
if jacobianAvg < tolerance.jacobianTolerance
    disp("Jacobian less than tolerance so we may have some singularities");
    % Regridd boolean
    regridBool = "True";
    return;
end


end