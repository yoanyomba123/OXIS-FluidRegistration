function Tinterp = performLinearInterpolation(Template, DisplacementVectorField, gridObject)
uVecX = DisplacementVectorField.x;
uVecY = DisplacementVectorField.y;

uVecXLen = length(uVecX);
uVecYLen = length(uVecY);

dx = gridObject.dx;
dy = gridObject.dy;
X = 0; Y = 0;
Tx = Template;
Ty = Template;

for i = 1: uVecXLen - 1 
    for j = 1: uVecYLen - 1
        X(i, j) = abs(gridObject.grid.x(i) + 1 - uVecX(i,j));
        Y(i, j) = abs(gridObject.grid.y(j) + 1 - uVecY(i,j));
        
       
    end
end

for i = 1: uVecXLen - 1
    for j = 1: uVecYLen -1
        newX = floor(X(i,j));
        newY = floor(Y(i,j));
        if newY == 0
            newY = 1;
        end
        
        if newX == 0
            newX = 1;
        end
        
        Tx(i,j) = Template(newX,j) + ((Template(newX+1,j)-Template(newX,j))/dx) * (X(i,j)-newX);
        Ty(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);
    end
end

% apply boundary conditions
Tinterp = struct();
Tinterp.x = applyDirichletBC(Tx);
Tinterp.y = applyDirichletBC(Ty);
end


function outputImage = applyDirichletBC(Template)
    [rows, cols] = size(Template);
    for i = 1: rows
        Template(i, 1) = 0;
        Template(i, end) = 0;
    end
    
    for j = 1: cols
       Template(1, j) = 0;
       Template(end , j) = 0;
    end
    
    outputImage = Template;
end
