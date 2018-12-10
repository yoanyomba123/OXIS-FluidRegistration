function outputImage = generateTrueGridImage(X,Y,Image)
% Obtain True Template
for i = 2:length(X)
    for j = 2:length(Y)
        outputImage(i,j) = Image(X(i), Y(j));
    end
end

end

