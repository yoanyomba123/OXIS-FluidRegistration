function [Tinterp, U] = performLinearInterpolation(Template, tRegrid,U, gridObject)
uVecX = U.x;
uVecY = U.y;

uVecXLen = length(uVecX);
uVecYLen = length(uVecY);

dx = gridObject.dx;
dy = gridObject.dy;
X = 0; Y = 0;
Tx = Template;
Ty = Template;
Tinterp = Template;

%U = applyDirichletBC(U);

for i = 1: uVecXLen - 1 
    for j = 1: uVecYLen - 1
        X(i, j) = gridObject.grid.x(i) - U.x(i,j);
        Y(i, j) = gridObject.grid.y(j) - U.y(i,j);
    end
end

for i = 1: uVecXLen - 1
    for j = 1: uVecYLen -1
        newX = floor(X(i,j));
        % interpolate in the x direction
        if(  newX > 0 & newX < length(Template) )
            Tx(i,j) = Template(newX,j) + ((Template(newX+1,j)-Template(newX,j))/dx) * (X(i,j)-newX);
        end
    end
end

for i = 1: uVecXLen - 1
    for j = 1: uVecYLen -1
        newY = floor(Y(i,j));
        % interpolate in the y direction
        if( newY > 0  & newY < length(Template))
            Ty(i,j) = Tx(i,newY) + ((Tx(i,newY+1)-Tx(i,newY))/dy) * (Y(i,j)-newY);
        end
    end
end


% apply boundary conditions
Tinterp = applyNeumanBC(Ty, gridObject.width);
end


function Image = applyNeumanBC(Image, width)
    [rows, cols] = size(Image);
    for j = 1:cols
        for i = 1:width
            Image(i,j) = 0;% Image(width,j);
        end
        
        for i=rows-width:rows
            Image(i,j) = 0;%Image((rows-1)-width,j);	
        end
    end
    
    for i=1:rows
		for j = 1:width
            Image(i,j) = 0;%Image(i,width);
        end
		for j = cols-width:cols
            Image(i,j) = 0;%Image(i,(cols-1)-width);
        end
    end
    
    for i = 1:width
		Image(  i  ,  i  ) = 0;% Image( width     ,   width    );
		Image(rows-1-i,  i  ) = 0;%Image((rows-1)-width,   width    );
		Image(  i  ,cols-1-i) = 0;%Image( width     , (cols-1)-width);
		Image(rows-1-i,cols-1-i) = 0;%Image((rows-1)-width, (cols-1)-width);
    end
   
    
%     for i = 1: rows
%         Image.x(i,1) = 0;
%         Image.x(i, end) = 0;
%     end
%     
%     for j = 1: cols
%        Image.y(1,j) = 0;
%        Image.y(end, j) = 0;
%     end
    
    
end

