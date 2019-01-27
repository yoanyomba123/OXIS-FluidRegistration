function MSE= computeMeanSquareError(Template, Source)
[M, N] = size(Template);
error = Template - (Source);
MSE = sum(sum(error .* error)) / (M * N);
disp(MSE);
end

