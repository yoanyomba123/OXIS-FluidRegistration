function [Tinterp, U] = performLinearInterpolation(Template, tRegrid,U, gridObject)
uVecX = tRegrid.x;
uVecY = tRegrid.y;

uVecXLen = length(uVecX);
uVecYLen = length(uVecY);

dx = gridObject.dx;
dy = gridObject.dy;
X = 0; Y = 0;
Tx = Template;
Ty = Template;

U = applyDirichletBC(U);

for i = 1: uVecXLen - 1 
    for j = 1: uVecYLen - 1
        X(i, j) = abs(gridObject.grid.x(i) + 1 + U.x(i,j));
        Y(i, j) = abs(gridObject.grid.y(j) + 1 + U.y(i,j));
    end
end

for i = 1: uVecXLen - 1
    for j = 1: uVecYLen -1
        newX = floor(X(i,j));
        newY = floor(Y(i,j));
        if newY == 0
            newY = 1;
        end
        
        if newY >= length(Template)
            newY = length(Template)-1;
        end
        
        if newX == 0
            newX = 1;
        end
        
        if newX >= length(Template)
            newX = length(Template) - 1;
        end
        
        Tx(i,j) = Template(newX,j) + ((Template(newX+1,j)-Template(newX,j))/dx) * (X(i,j)-newX);
        Ty(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);
    end
end

% apply boundary conditions
Tinterp = struct();
Tinterp.x = Tx;
Tinterp.y = Ty;

end


function displacementField = applyDirichletBC(U)
    [rows, cols] = size(U.x);
    for i = 1: rows
        U.x(i,1) = 0;
        U.x(i, end) = 0;
    end
    
    for j = 1: cols
       U.y(1,j) = 0;
       U.y(end, j) = 0;
    end
    
    displacementField = U;
end
