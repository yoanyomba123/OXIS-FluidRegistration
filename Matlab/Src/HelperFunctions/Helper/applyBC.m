function displacementField = applyBC(U)
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
