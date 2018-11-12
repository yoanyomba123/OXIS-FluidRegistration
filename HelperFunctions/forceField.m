function force = forceField(Template, Source, x, y, Ux, Uy, dx, dy);
% Template and Source are images
% x is a vector(point) comprised of two values delineating
% the x component and the y component of the point
% u is a direction (dispalcement) vector comprised of
% two components (x and y)

% Computes Force Field
SSD = 0.0;
alpha = 0.5;

% Want to start iteration at points not on edges
xLen = length(x);
yLen = length(y);
% Obtain the body field
for i=2:xLen-1;
    for j=2:yLen-1;
        xPos = [x(i), y(j)]; % obtain position
        xComp = xPos(1);
        yComp = xPos(2);
        TemplatePos = double(Template(xComp-Ux(xComp), yComp-Uy(yComp)));
        SourcePos = double(Source(xComp,yComp));
        SSD = SSD + (alpha/2)*((TemplatePos - SourcePos) ...
        *(TemplatePos - SourcePos));
    end
end

% scale the Sum of Squared Difference
SSD = SSD/((xLen-2 + yLen-2));

% Now Compute Force Field
for i=2:xLen;
    for j=2:yLen;
        xPos = [x(i), y(j)]; % obtain position
        xComp = xPos(1);
        yComp = xPos(2);
        TemplateXComp = xComp - Ux(xComp);
        TemplateYComp = yComp - Uy(yComp);
        % compute the partial of T with respect to x and Y by Forward Difference
        % Scheme
        dTdx = (double(Template(TemplateXComp+1,TemplateYComp)) - double(Template(TemplateXComp-1,TemplateYComp))) / dx; % Partial of T with respect to x 
        dTdy = (double(Template(TemplateXComp,TemplateYComp+1)) - double(Template(TemplateXComp,TemplateYComp-1))) / dy; % Partial of T with respect to y
        % Compute the actual force field by taking the difference between
        % the template image and the source multiplied by the computed
        % gradient
        fx(i,j) = -alpha*(double(Template(TemplateXComp,TemplateYComp)) - double(Source(i,yComp))) * dTdx;
		fy(i,j) = -alpha*(double(Template(TemplateXComp,TemplateYComp)) - double(Source(i,yComp))) * dTdy;
    end
end

fx = (-1.0).*fx;
fy = (-1.0).*fy;

force(:,:,1) = fx;
force(:,:,2) = fy;
end
