function regridBool = computeJacobian(gridObject, U, maxPertubation, perturbation, tolerance)
regridBool = "false";
dx = gridObject.dx;
dy = gridObject.dy;

x = gridObject.x;
y = gridObject.y;

dx2 = dx * 2.0;
dy2 = dy * 2.0;

% expX =  U.x - (maxPertubation .* perturbation.x);
% expY =  U.y - (maxPertubation .* perturbation.y);
% Jac = zeros(length(U.x), length(U.y));
% 
% for i = 2: length(expX)-1
%     for j = 2: length(expY)-1
%         dex_dx = (expX(i+1,j) - expX(i-1, j))/2;
%         dex_dy = (expX(i, j+1) - expX(i, j-1))/2;
%         
%         dey_dx = (expY(i+1,j) - expY(i-1, j))/2;
%         dey_dy = (expY(i, j+1) - expY(i, j-1))/2;
%         
%         Jac(i, j) = ((dey_dy .* dex_dx) - (dey_dx .* dex_dy));
%     end
% end
% process the jacobian
%jacobianObject = struct();
%jacobianObject.referenceFieldEntityX = U.x; %gridObject.grid.x - U.x - (maxPertubation .* perturbation.x);
%jacobianObject.referenceFieldEntityY = U.y; %gridObject.grid.y - U.y - (maxPertubation .* perturbation.y);

%entityX = jacobianObject.referenceFieldEntityX;
%entityY = jacobianObject.referenceFieldEntityY;

%entityX = padarray(entityX, [1,1]); entityY = padarray(entityY, [1,1]);
%jacobianAvg = 0.0;


for i = 2:length(U.x)-2
    for j= 2:length(U.y)-2           
        % process the jacobian
        dux_dx = (U.x(i+1, j) - U.x(i-1, j))/2;
        dux_dy = (U.x(i, j+1) - U.x(i, j-1))/2;

        duy_dy = (U.y(i, j+1) - U.y(i, j-1))/2;
        duy_dx = (U.y(i+1, j) - U.y(i-1, j))/2;
        
        rx_dx = maxPertubation * (perturbation.x(i+1, j) - perturbation.x(i-1, j));
        rx_dy = maxPertubation * (perturbation.x(i, j+1) - perturbation.x(i, j-1));
        ry_dy = maxPertubation * (perturbation.y(i, j+1) - perturbation.y(i, j+1));
        ry_dx = maxPertubation * (perturbation.y(i+1, j) - perturbation.y(i-1, j));

        
        a = dx - dux_dx - rx_dx;
        b = dx - duy_dx - ry_dx;
        c = dy - dux_dy - rx_dy;
        d = dy - duy_dy - ry_dy;
        
        
        %Jacobian = [dx_dx, dx_dy; dy_dx, dy_dy];
        %jacobianAvg = jacobianAvg + abs(det(Jacobian));
        
        % store jacobian determinant
        Jac(i, j) = ((a .* d) - (b * c));
    end
end

if abs(min(min(Jac))) < tolerance.jacobianTolerance
    disp(abs(min(min(Jac))));
    disp("Jacobian less than tolerance so we may have some singularities");
    % Regridd boolean
    regridBool = "True";
    return;
end


end

