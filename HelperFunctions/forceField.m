function force = forceField(Template, Source, displacement, gridObject, method);
% Template and Source are images
% x is a vector(point) comprised of two values delineating
% the x component and the y component of the point
% u is a direction (dispalcement) vector comprised of
% two components (x and y)

% Computes Force Field
SSD = 0.0;
alpha = 0.5;

x = gridObject.x;
y = gridObject.y;
Ux = displacement.x;
Uy = displacement.y;
dx = gridObject.dx;
dy = gridObject.dy;

% Want to start iteration at points not on edges
xLen = length(x);
yLen = length(y);
% % Obtain the body field
% for i=2:xLen-1;
%     i = real(i);
%     for j=2:yLen-1;
%         j = real(j);
%         xPos = [x(i), y(j)]; % obtain position
%         xComp = ceil(xPos(1));
%         yComp = ceil(xPos(2));
%         TemplatePos = double(Template(ceil(xComp-Ux(xComp)), ceil(yComp-Uy(yComp))));
%         SourcePos = double(Source(xComp,yComp));
%         SSD = SSD + (alpha/2)*((TemplatePos - SourcePos) ...
%         *(TemplatePos - SourcePos));
%     end
% end
% 
% % scale the Sum of Squared Difference
% SSD = SSD/((xLen-2 + yLen-2));

% Now Compute Force Field
for i=2:xLen;
    for j=2:yLen;
        xPos = [x(i), y(j)]; % obtain position
        xComp = xPos(1);
        yComp = xPos(2);
        TemplateXComp = ceil(real(ceil(xComp - Ux(xComp))));
        TemplateYComp = ceil(real(ceil(yComp - Uy(yComp))));
        
        % compute the partial of T with respect to x and Y by Forward Difference
        % Scheme
        dTdx = (double(Template.x(TemplateXComp+1,TemplateYComp)) - double(Template.x(TemplateXComp-1,TemplateYComp))) / dx; % Partial of T with respect to x 
        dTdy = (double(Template.y(TemplateXComp,TemplateYComp+1)) - double(Template.y(TemplateXComp,TemplateYComp-1))) / dy; % Partial of T with respect to y
        % Compute the actual force field by taking the difference between
        % the template image and the source multiplied by the computed
        % gradient
        fx(i,j) = -alpha*(double(Template.x(TemplateXComp,TemplateYComp)) - double(Source(TemplateXComp,TemplateYComp))) * dTdx;
		fy(i,j) = -alpha*(double(Template.y(TemplateXComp,TemplateYComp)) - double(Source(TemplateXComp,TemplateYComp))) * dTdy;
        
        % TODO : FIX THIS LATER
        Sample = Template.x - Source;
        
        
        if(method == "QN")
           if i == j
               Data = Sample(1:i, 1:j); 
               centered_diff = full(gallery('tridiag', i, 1,-2,1));
               Hessian = Data .* centered_diff .* centered_diff;
               Gradient = gradient(Data);
               if(det(Hessian) > 0)
                fx(1:i, 1:j) = fx(1:i, 1:j) - (inv(Hessian) .* Gradient);
                fy(1:i, 1:j) = fy(1:i, 1:j) - (inv(Hessian) .* Gradient);
               end
           end
        end
    
    end
end

%fx = (-1.0).*fx;
%fy = (-1.0).*fy;
force = struct();
force.x = fx;
force.y = fy;
end
