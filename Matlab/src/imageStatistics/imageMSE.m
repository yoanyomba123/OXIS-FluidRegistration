function mse = imageMSE(image1,image2)

if (size(image1) ~= size(image2))
    disp("Images Must Be Of Same Dimensions");
    return
end

[rows,cols] = size(image1);

sum = 0.0;

for i = 1: length(image1)
    for j = 1: length(image2)
        diff = image1(i,j) - image2(i,j);
        sum = sum + (diff * diff);
    end
end

mse = sum / (rows * cols);
end

