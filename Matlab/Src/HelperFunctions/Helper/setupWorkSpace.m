% Method initializes various objects to be used in registration and stores
% them in a .mat file to be loaded into calling method workspace
function setupWorkSpace(Template, Source, numPoints, iter, mu, lambda, deformationTol, jacobianTol, distanceTol)
    [rows, cols] = size(Template);
    
    Template = Template;
    Source = Source;
    
    maxIter = iter;
    params = struct();
    params.mu = mu;
    params.lambda = lambda;
    
     % define stencils
    stencil = struct();
    stencil.S11 = [0, (params.lambda + 2*params.mu), 0;
           params.mu, -2*(params.lambda + 3*params.mu), params.mu;
           0, (params.lambda+2*params.mu), 0 ];
    stencil.S12 = [1, 0, -1;
           0, 0 ,0;
           -1, 0, 1] .* ((params.lambda + params.mu)/4);
    stencil.S22 = stencil.S11';
    stencil.S21 = stencil.S12';
    
    % tolerance definition
    tolerance = struct();
    tolerance.deformationTolerance = deformationTol;
    tolerance.jacobianTolerance = jacobianTol;
    tolerance.distanceTolerance = distanceTol;
    tolerance.regridTolerance = 20;
    tolerance.mse = 1e-13;
    
    gridObject = struct();
    gridObject.numXPoints = numPoints;
    gridObject.numYPoints = numPoints;
    gridObject.grid = struct();

    % generate points that are not on the boundary of the image
    x = linspace(0, rows-2, gridObject.numXPoints); x = ceil(x);
    y = linspace(0, cols-2, gridObject.numYPoints); y = ceil(y);
    [X,Y] = meshgrid(x, y);
    X = X + 1; Y = Y + 1;
    x = x + 1; y = y + 1;
    gridObject.rows = rows;
    gridObject.cols = cols;
    gridObject.x = x;
    gridObject.y = y;
    gridObject.grid.x = X;
    gridObject.grid.y = Y;
    gridObject.width = 1;
    gridObject.dx = 1*ceil(gridObject.rows/gridObject.numXPoints);
    gridObject.dy = 1*ceil(gridObject.cols/gridObject.numYPoints);
    
    U = struct();
    U.x = zeros(gridObject.numXPoints, gridObject.numYPoints);
    U.y = zeros(gridObject.numXPoints, gridObject.numYPoints);
    
     % initialze regrid components
    structVals = {'x', 'y', 'template'};
    yQ = cell(1, maxIter);
    yRegrid = struct();
    yRegrid.x = zeros(gridObject.numXPoints, gridObject.numYPoints);
    yRegrid.y = zeros(gridObject.numXPoints, gridObject.numYPoints);
    yRegrid.template = [];
    regridCounter = 1;
    yQ{regridCounter} = yRegrid;

    tK = cell(1, maxIter);
    tRegrid = struct();
    tRegrid.x = zeros(gridObject.numXPoints, gridObject.numYPoints);
    tRegrid.y = zeros(gridObject.numXPoints, gridObject.numYPoints);

    wK = cell(1, maxIter);
    wRegrid = struct();
    wRegrid.x = zeros(gridObject.numXPoints, gridObject.numYPoints);
    wRegrid.y = zeros(gridObject.numXPoints, gridObject.numYPoints);
    
    gridObject.sampleTemplate = generateTrueGridImage(gridObject.x, gridObject.y, Template);
    gridObject.sampleSource = generateTrueGridImage(gridObject.x, gridObject.y, Source);
    save("variables");
end
