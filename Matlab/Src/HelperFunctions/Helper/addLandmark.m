function outputImage = addLandmark(inputImage, numLandMarks)

% obtain image szie
[rows, cols] = size(inputImage);

% obtain range of operation on image
xmin = 1;
xmax = cols-2;
ymin = 1;
ymax = rows-2;

% set landmarks
for i = 1:numLandMarks
    x = xmin + rand(1, 1) * (xmax - xmin);
    y = ymin + rand(1,1) * (xmax - xmin);
    
    x = ceil(x);
    y = ceil(y);
    inputImage = insertMarker(inputImage, [x,y], 'o', 'color', 'black', 'size', 10);
end

outputImage = mat2gray(inputImage);

end

