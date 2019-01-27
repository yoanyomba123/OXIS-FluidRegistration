function [Template, Source, preMSE, postMSE, preIntDiff, postIntDiff, TemplateOut, gridObject, U] = fluidRegistration2DImpl(Template, Source, numpoints, iter, mu, lambda, tolerance)

% setup workspace environment
setupWorkSpace(Template, Source, numpoints, iter, mu, lambda, tolerance.deformationTolerance, tolerance.jacobianTolerance, tolerance.distanceTolerance);

% load workspace variables
load("variables.mat");

i = 1;

% Obtain Initial Interpolated Template
[interpT, U] = performLinearInterpolation(Template,Source,U,gridObject);

% Define The Jacobian Matrix
Jacobian = zeros(length(x), length(y));

numRegrids = 0;

% dmin = immse(Template, Source);

% compute the mse of the two images prior to registration
preMSE = immse(Template, Source);

% compute the intensity diff of the two images prior to registration
preIntDiff = imabsdiff(Template, Source);

initialMSE = 1;

while i < iter
    % store prior tolerance value for optimization
    deformationDistTolprevious = tolerance.deformationTolerance;
    
    % perform 2D interpolation
    % turn this into a function
%     wRegrid.x = interpn(yQ{regridCounter}.x,gridObject.grid.x - U.x,'linear');
%     wRegrid.y =  interpn(yQ{regridCounter}.y,gridObject.grid.y - U.y,'linear');
%     wK{i} = wRegrid;
    
    tRegrid.x = X - U.x;
    tRegrid.y = Y - U.y;
    [tK{i}, U] = performLinearInterpolation(Template,tRegrid,U,gridObject);

    % Minimization is performed in the forcefield function
    force = computeForceFieldSSD(tK{i}, Source, U, gridObject);
    
    % evaluate the velocity vector fields by solving the linear discretized
    % PDE
    V = computeVelocityVectorFields(stencil, force, gridObject);
    
    % compute the pertubation of the displacement field
    [perturbation, delta] = computePertubation(gridObject, U, V, tolerance);    
    
    % compute the jacobian and update u as well as perform regridding if
    % necessary
    regridBool = computeJacobian(gridObject, U, delta, perturbation, tolerance);
    
    if(numRegrids <= tolerance.regridTolerance)
        if(regridBool == "True" && i > 1)
            count = 0;
            while(regridBool == "True")
               numRegrids = numRegrids + 1;
               % lower step size
               Template = deformTemplate(Template, U, gridObject);
               U = clearDisplacement(U);
               tolerance.distanceTolerance = tolerance.distanceTolerance * 2;
               i = i + 1;
               [perturbation, delta] = computePertubation(gridObject, U, V, tolerance); 
               regridBool = computeJacobian(gridObject, U, delta, perturbation, tolerance);
               count = count + 1;
               
               if(numRegrids > tolerance.regridTolerance)
                  break; 
               end
            end
            % deform template
           Template = deformTemplate(Template, U, gridObject);
           
           figure; imagesc(Template);
           
           U = clearDisplacement(U);
        else
            % apply pertubation to displacement field
            if(delta > 0)
               % perturb the displacement field 
               U.x = U.x + (delta .* perturbation.x); 
               U.y = U.y + (delta .* perturbation.y);
            end
            
            % obtain difference between source and template
            diff = immse(deformTemplate(Template, U, gridObject), Source);
%             if(diff < dmin)
%                 % increase the tolerance per iteration
%                 tolerance.distanceTolerance = tolerance.distanceTolerance * 2;
%             else
%                 tolerance.distanceTolerance = tolerance.distanceTolerance/2;
%             end
            
        end
     else
        disp("Max Amount Of Regrids Encountered")
        break;   
    end
    
%     currentMSE = immse(Template, Source);
%     if(i > 1)
%        if norm((Source - Template) - (Source  - Template),2) <=  norm((Source  - Template) .* tolerance.distanceTolerance,2)
%             return;
%        end
%     end
%     
%      % terminating condition (Max iteration reached or mse tolerance 
%     mseRateOfChange = abs((currentMSE - initialMSE)/(initialMSE))
%     if i > maxIter | mseRateOfChange <= tolerance.mse
%        return
%     end
    
    if i > iter
       return 
    end
    
%     initialMSE = currentMSE;
    i = i + 1;
end

% Warp The Image
TemplateOut = 0;
for d = 1: length(U.x)
    for j = 1: length(U.y)
        if(x(d) - U.x(d, j)) <= length(Template) & (y(j) - U.y(d, j)) <= length(Template) & (x(d) -  U.x(d, j)) > 0 & (y(j) - U.y(d, j)) > 0
            TemplateOut(ceil(x(d) -  U.x(d, j)),ceil(y(j) -  U.y(d, j))) =  Source(x(d),y(j)); 
        end
    end
end

gridObject.grid.x = gridObject.grid.x - U.x;
gridObject.grid.y = gridObject.grid.y - U.y;

postMSE = immse(Source, TemplateOut);
postIntDiff = imabsdiff(Source, TemplateOut);
end

