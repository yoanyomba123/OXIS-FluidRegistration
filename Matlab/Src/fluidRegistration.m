% Function commences registration of images by use of the fluid registration algorithm
% @param - path to a directory containing dicom images
% @result - a structure containing all registration results 
%         - structure contains the Template and Source Images, pre and post Registration pixelwise MSE
%           the Deformed image, x and y coordinate planes, the displacement fields, and a gridObject

function regRes = fluidRegistration(dir)
    regRes = struct();
    files = readImagesFromDirectory(dir);

    
    iter = 500;
    mu = 10; lambda = 10;

    % tolerance definition
    tolerance = struct();
    tolerance.deformationTolerance = 50;
    tolerance.jacobianTolerance = 0.5;
    tolerance.distanceTolerance = 1e-2;
    tolerance.mse = 1e-12;

    % read all images and augment them then register
    for i = 1:2;
        numpoints = 510;
        % transform / warp image to create template
        Template = imrotate(files{i},-60,'bilinear','crop'); % rotate the template image

        % obtain source image
        Source = files{i};
        
        % perform registration of images
        [Template, Source, preMSE, postMSE, preIntDiff, postIntDiff, TemplateOut, gridObject, U] = fluidRegistration2DImpl(Template, Source, numpoints, iter, mu, lambda, tolerance);

        % store relevant stats in a struct for each image
        % mse pre and post registration, intensity difference pre and post
        % registration, deformed image, as well as template and source images
        regRes(i).Template = Template;
        regRes(i).Source = Source;
        regRes(i).preRegMSE = preMSE;
        regRes(i).postRegMSE = postMSE;
        regRes(i).preIntDiff = preIntDiff;
        regRes(i).postIntDiff = postIntDiff;
        regRes(i).deformedImage = TemplateOut;
        regRes(i).x = gridObject.grid.x;
        regRes(i).y = gridObject.grid.y;
        regRes(i).gridData = gridObject;
        regRes(i).displacement = U;
        regRes(i).preCorrCoeff = corr2(Template, Source);
        regRes(i).postCorrCoeff = corr2(TemplateOut, Source);
    end
end
