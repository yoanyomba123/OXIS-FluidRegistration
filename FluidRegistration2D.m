%% Viscous Fluid Image Registrion
% Author: D Yoan L. Mekontchou Yomba
% Date: 11/16/2018
% Purpose:  
%   The purpose of this script is to implement a fluid registration model
%   in matlab. The Registration model makes use of the Eularian reference
%   frame and uses the sum of square difference as a cost function.
%% Clear Up Workspace
clc; clear; close all;

%% Add all paths to current workspace recursively
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Load In Images Into Workspace
[Template, Source] = loadImages("Data");

Template =  Template;
Source =  Source;


Template = imrotate(Template,-5,'bilinear','crop'); % rotate the template image

% display and visualize both images
figure; imagesc([Template, Source]);

% display and visualize difference between both images
 figure; imshowpair(Template, Source);
maxdiff = max(max(Template - Source))
%% Define Initial Conditions 

% params
params = struct();
params.mu = 1;
params.lambda = 1;

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
tolerance.deformationTolerance = 300;
tolerance.jacobianTolerance = 0.35;
tolerance.distanceTolerance = 1e-7;
tolerance.mse = 1e-10;

% max iteration terminating condition
maxIter = 100;

% grid definition
[rows, cols] = size(Template);
gridObject = struct();
gridObject.numXPoints = 50;
gridObject.numYPoints = 50;
gridObject.grid = struct();

% generate points that are not on the boundary of the image
x = linspace(0, rows-2, gridObject.numXPoints); x = ceil(x);
y = linspace(0, cols-2, gridObject.numYPoints); y = ceil(y);
[X,Y] = meshgrid(x, y);
gridObject.rows = rows;
gridObject.cols = cols;
gridObject.x = x;
gridObject.y = y;
gridObject.grid.x = X;
gridObject.grid.y = Y;
gridObject.width = 1;
gridObject.dx = 1*ceil(gridObject.rows/gridObject.numXPoints);
gridObject.dy = 1*ceil(gridObject.cols/gridObject.numYPoints);

% initilaze displacement field
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

% central Diffence Matrix Operator
centralDiffMatOperator = full(gallery('tridiag', length(gridObject.x), -1,2, -1));

% define fourier matrix and inverse fourier matrix operator
fftMatOperator = dftmtx(gridObject.numXPoints);
fftMatInvOperator = inv(fftMatOperator);

A = struct();
V = struct();
Jacobian = struct();
deltalU = struct();

% Obtain Template and Source Based On Grid Point Definitions 
gridObject.sampleTemplate = generateTrueGridImage(gridObject.x, gridObject.y, Template);
gridObject.sampleSource = generateTrueGridImage(gridObject.x, gridObject.y, Source);


wk = struct();
displacementVector = cell(1, maxIter);
x = x + 1;
y = y + 1;
TemplateSet = cell(1, maxIter);
%%
figure;
priorDelta = 0;
i = 1;
initialMSE = 1;
while 1;
    % store prior tolerance value for optimization
    deformationDistTolprevious = tolerance.deformationTolerance;
    
    % perform 2D interpolation
    % turn this into a function
    wRegrid.x = interpn(yQ{regridCounter}.x,gridObject.grid.x - U.x,'linear');
    wRegrid.y =  interpn(yQ{regridCounter}.y,gridObject.grid.y - U.y,'linear');
    wK{i} = wRegrid;
    
    tRegrid.x = X - wRegrid.x - U.x;
    tRegrid.y = Y - wRegrid.y - U.y;
    [tK{i}, U] = performLinearInterpolation(Template,tRegrid,U,gridObject);

    % Minimization is performed in the forcefield function
    force = forceField(tK{i}, Source, U, gridObject, "none");
    
    drawnow
    % visualize the force field on the image
    %visualize(force.x, force.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    %disp("Displaying force fields");
    %pause(2);
    
    % obtain V
    % Turn this into a function
    Sx = stencil.S11 + stencil.S12;
    Sy = stencil.S21 + stencil.S22;
    
    A.x = conv2(gridObject.sampleTemplate, Sx, 'same');
    A.y = conv2(gridObject.sampleTemplate, Sy, 'same');
    A.FFTx = real(fft2(A.x));%real(fftMatOperator .* A.x .* fftMatInvOperator);
    A.FFTy = real(fft2(A.y)); %real(fftMatOperator .* A.y .* fftMatInvOperator);
    
    Dx = pinv(A.FFTx);
    Dy = pinv(A.FFTy);
    
    V.x = (Dx) .* force.x;
    V.y = (Dy) .* force.y;
    V = applyBC(V);

    % visualize the velocity field on the image
    %visualize(V.x, V.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    %disp("Displaying Velocity Vector fields");
    %pause(2);
 
   [pertubation, delta] = computePertubationAndUpdateDisplacement(gridObject, U, V, tolerance, regridCounter, wK, yQ, tK, i);
   
   regridBool = computeJacobian(gridObject, U, delta, pertubation, tolerance);
   if regridBool == "false"
        U.x = U.x + (pertubation.x .* delta);
        U.y = U.y + (pertubation.y .* delta);
   else
        singularityCount = 0;
        while regridBool == "True"
           [regridCounter, wK, U, yQ, tK, pertubation, delta, iteration] = perfromRegridding(gridObject, regridBool, regridCounter, wK, U, yQ, tK, pertubation, delta, i, Template);
           [pertubation, delta] = computePertubationAndUpdateDisplacement(gridObject, U, V, tolerance, regridCounter, wK, yQ, tK, i);
           regridBool = computeJacobian(gridObject, U, delta, pertubation, tolerance);
           singularityCount = singularityCount + 1;
           
           if(singularityCount > 5)
               break;
           end
        end
    end
  
    %visualize(U.x, U.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    %disp("Displaying Displacement Vector fields");
    %pause(2);
    
    currentMSE = immse(tK{i}, Source);
    if(i > 1)
       if norm((Source - tK{i-1}) - (Source  - tK{i}),2) <=  norm((Source  - tK{i}) .* tolerance.distanceTolerance,2)
            return;
       end
    end
    
     % terminating condition (Max iteration reached or mse tolerance 
    mseRateOfChange = abs((currentMSE - initialMSE)/(initialMSE))
    if i > maxIter | mseRateOfChange <= tolerance.mse
       return
    end
    
    initialMSE = currentMSE;
    i = i +1;
end
%%
TemplateOut = 0;
SourceOut = Source;
% U.x = U.x ./ tolerance.distanceTolerance;
% U.y = U.y ./ tolerance.distanceTolerance;

c = 0;
for d = 1: length(U.x)
    for j = 1: length(U.y)
        if(x(d) - U.x(d, j)) <= length(Template) & (y(j) - U.y(d, j)) <= length(Template) & (x(d) - U.x(d, j)) > 0 & (y(j) - U.y(d, j)) > 0
            TemplateOut(ceil(abs(x(d) -  U.x(d, j))), ceil(abs(y(j) -   U.y(d, j)))) = Template(x(d),y(j)); 
            SourceOut(x(d), y(j)) = Source(x(d), y(j));
            c= c + 1;
        end
    end
end

figure; imagesc([TemplateOut]); colormap gray;
title("Deformed Template (Iteration: 85)");


% figure; imshowpair(TemplateOut, Source);
% title("Deformed Template Vs Source (Iteration: 85)");
% 
figure; imagesc([TemplateOut Template]);
% title("Deformed Template vs Original Template (Iteration: 85)")
%%
U.x =  (wK{i}.x + real((U.x)));
U.y =  (wK{i}.y + real((U.y)));

field(:,:,1) = ceil(X - U.x);
field(:,:,2) = ceil(Y - U.y);

%Tout = imwarp(Template, field);
x = x+1;
y = y+1;
templateT = 0;
templateInit = 0;
SourceGrid = 0;
for i = 1: length(U.x)
    for j = 1: length(U.y)
        TemplateT(ceil(abs(x(i) + U.x(i, j))), ceil(abs(y(j) + U.y(i, j)))) = Template(x(i),y(j)); 
        TemplateInit(x(i),y(j)) = Template(x(i),y(j));
    end
end
tout = imwarp(Template, field);


figure; imshowpair(TemplateInit, SourceGrid,'ColorChannels', 'green-magenta');
figure; imshowpair(TemplateT, SourceGrid, 'ColorChannels', 'red-cyan');

%%
displacementfield(:,:,1) = ceil(X + U.x);
displacementfield(:,:,2) = ceil(Y + U.y);

figure;plot(displacementfield(:,:,1),displacementfield(:,:,2)); hold on;

visualize(U.x, U.y, displacementfield(:,:,1), displacementfield(:,:,2), gridObject.sampleTemplate);
disp("Displaying Displacement Vector fields");
pause(1);

%%
output = imwarp(sampleTemplate, displacementfield);
figure; imshowpair(output,sampleTemplate);
title("Registered Image vs Template");

figure; imshowpair(output, sampleSource,"ColorChannels", 'red-cyan');
title("difference between the Registered Image amd Source");

figure; imshowpair(output, sampleTemplate,"ColorChannels", 'red-cyan');
title("difference between the Output Image amd Template");


RegisteredImage = zeros(length(U.x), length(U.y));

for i = 1: length(U.x)
    for j =1: length(U.y)     
        x_hat = U.x(i,j);
        y_hat = U.y(i,j);
        x_def = ceil(i + x_hat);
        y_def = ceil(j + y_hat);
        if(x_def > 0 && y_def > 0)
           RegisteredImage(x_def, y_def) = sampleTemplate(i,j);
        end
    end
end

RegisteredImage = RegisteredImage(1:end-1,1:end-1);
figure; imagesc(RegisteredImage);
figure; imagesc(sampleTemplate);



% figure; imagesc([RegisteredImage sampleTemplate]);
% figure; imshowpair(RegisteredImage, sampleSource,'ColorChannels','red-cyan');
figure; imshowpair(RegisteredImage, sampleTemplate,'ColorChannels','red-cyan');

% tOut(:,:,1) = interp2(Template,  gridObject.grid.x - U.x);
% tOut(:,:,2) = interp2(Template, gridObject.grid.y - U.y);
% 
% for i = 2: length(U.x)
%     for j = 2: length(U.y);
%         output(i, j) = Template(x(i), y(j)) + U.x(i, j)+ U.y(i, j) ;
%     end
% end
