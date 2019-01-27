function deformedImage = deformImage(Template, DisplacementField, X, Y)
len = length(X);

deformedImage = Template;
for i = 1: len
    for j = 1: len
        newX = (X(i)-DisplacementField.x(i,j));
        newY = (Y(i)-DisplacementField.y(i,j));
        if newX > 0 & newY > 0 & newX <= rows & newY <= cols
            deformedImage(newX,newY) = deformedImage(X(i), Y(j));
        end
    end
end

end

