function outputImage = performLinearInterpolation(Template, DisplacementVectorField, gridObject)
uVecX = DisplacementVectorField.x;
uVecY = DisplacementVectorField.y;

uVecXLen = length(uVecX);
uVecYLen = length(uVecY);

dx = gridObject.dx;
dy = gridObject.dy;
X = 0; Y = 0;
Tx = Template;
Ty = Template;
Tinterp = Template;

for i = 1: uVecXLen - 1 
    for j = 1: uVecYLen - 1
        X(i, j) = i - uVecXlen(i,j);
        Y(i, j) = j - uVecYlen(i,j);
    end
end

for i = 1: uVecXLen - 1
    for j = 1: uVecYLen -1
        newX = long(floor(X(i,j))); newY = long(floor(Y(i,j)));
        Tx(i,j) = Template(newX,j) + ((Template(newX+1,j)-Template(newX,j))/dx) * (X(i,j)-newX);
        Ty(i,j,k) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);
    end
end

% apply boundary conditions

end

