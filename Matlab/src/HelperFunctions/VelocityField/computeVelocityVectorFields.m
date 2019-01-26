function V = computeVelocityVectorFields(stencil, force, gridObject)
    
    % Perform Stencils Mathematics with Stencils
    Sx = stencil.S11 + stencil.S12;
    Sy = stencil.S21 + stencil.S22;
    
    % Convolve image with stencils
    A.x = conv2(gridObject.sampleTemplate, Sx, 'same');
    A.y = conv2(gridObject.sampleTemplate, Sy, 'same');
    
    % Apply fft transform
    A.FFTx = real(fft2(A.x));
    A.FFTy = real(fft2(A.y)); 
    
    % invert the circular matrix by use of perose-pseudoinverse
    Dx = pinv(A.FFTx);
    Dy = pinv(A.FFTy);
    
    % compute the velocity vector components
    V.x = (Dx) .* force.x;
    V.y = (Dy) .* force.y;
    V = applyBC(V);

end

