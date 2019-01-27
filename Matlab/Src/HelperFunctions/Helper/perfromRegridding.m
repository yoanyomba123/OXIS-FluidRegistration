function [regridCounter, wK, U, yQ, tK, pertubation, maxPertubation, iteration, Template] = perfromRegridding(gridObject, regridBool, regridCounter, wK, U, yQ, tK, pertubation, maxPertubation, iteration, Template)

% Perform Regriddind Conditionally
if(regridBool == "True")
    regridCounter = regridCounter+1;
    regridEntity = struct();
    regridEntity.x = wK{iteration}.x + U.x;
    regridEntity.y = wK{iteration}.y + U.y;
    yQ{regridCounter} = regridEntity;
    
    displacementfield(:,:,1) = gridObject.grid.x - U.x;
    displacementfield(:,:,2) = gridObject.grid.y - U.y;
    Template = imwarp(Template, displacementfield);
    U.x = U.x .* 0;
    U.y = U.y .* 0;
end

end

