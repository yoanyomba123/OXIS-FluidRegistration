%% Viscous Fluid Image Registrion
% Author: D Yoan L. Mekontchou Yomba
% Date: 11/16/2018
% Purpose:  
%   The purpose of this script is to implement a fluid registration model
%   in matlab. The Registration model makes use of the Eularian reference
%   frame and uses the sum of square difference as a cost function.

% Clear Up Workspace
clc; clear; close all;

% Add all paths to current workspace recursively
currentFolder = pwd;
addpath(genpath(currentFolder));

% Load In Images Into Workspace
[Template, Source] = loadImages("Data");

Template =  Template;
Source =  Source;

Template = imrotate(Template,-15,'bilinear','crop');
% Template  =  imrotate(Template,-10,'bilinear','crop');
% templateSourceDiff = Template - Source;

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
tolerance.deformationTolerance = 10;
tolerance.jacobianTolerance = 0.25;
tolerance.distanceTolerance = 100;
tolerance.delta = 0.00005;

maxIter = 100;

% grid definition
[rows, cols] = size(Template);
gridObject = struct();
gridObject.numXPoints = 40;
gridObject.numYPoints = 40;
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
% Obtain True Template
for i = 2:length(x)
    for j = 2:length(y)
        sampleTemplate(i,j) = Template(x(i), y(j));
        sampleSource(i,j) = Source(x(i), y(j));
    end
end
gridObject.sampleTemplate = sampleTemplate;
gridObject.sampleSource = sampleSource;
figure; imagesc([Template, Source]);
title("Template (Rotated) | Source")
% figure; imagesc([sampleSource]);
% title("Source Image");
% figure; imagesc([sampleTemplate]);
% title("Template Image");
% Algorithm
wk = struct();
displacementVector = cell(1, maxIter);
x = x + 1;
y = y + 1;
%%
figure;
priorDelta = 0;
i = 1;
while 1;
    % store prior tolerance value for optimization
    deformationDistTolprevious = tolerance.deformationTolerance;
    
    % perform 2D interpolation
    % turn this into a function
    wRegrid.x = interp2(yQ{regridCounter}.x,gridObject.grid.x + U.x,'linear');
    wRegrid.y =  interp2(yQ{regridCounter}.y,gridObject.grid.y + U.y,'linear');
    wK{i} = wRegrid;
    
    tRegrid.x = gridObject.grid.x + wRegrid.x + U.x;
    tRegrid.y = gridObject.grid.y + wRegrid.y + U.y;
    tK{i} = performLinearInterpolation(Template,tRegrid, gridObject);
%     for k = 1 : length(U.x)
%         for l = 1 : length(U.y)
%             tK{i}.x(x(k), y(l)) = tK{i}.x(x(k), y(l)) + U.x(k,l);
%             tK{i}.y(x(k), y(l)) = tK{i}.y(x(k), y(l)) + U.y(k,l);
% 
%         end
%     end

%     if(i > 1)
%        if max(max((Source - tK{i-1}.x) - (Source - tK{i}.y))) <=  max(max((Source - tK{i}.x).*tolerance.distanceTolerance)) | ...
%             max(max((Source - tK{i-1}.y) - (Source - tK{i}.y))) <=  max(max((Source - tK{i}.y).*tolerance.distanceTolerance))            
%             return;
%        end
%     end
    
    % Minimization is performed in the forcefield function
    force = forceField(tK{i}, Source, U, gridObject, "none");
    
%     drawnow
%     % visualize the force field on the image
    visualize(force.x, force.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    disp("Displaying force fields");
    pause(2);
    
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
    
%     % visualize the velocity field on the image
    visualize(V.x, V.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    disp("Displaying Velocity Vector fields");
    pause(2);
 
   [pertubation, delta] = computePertubationAndUpdateDisplacement(gridObject, U, V, tolerance, regridCounter, wK, yQ, tK, i);
   
   % terminating condition
    if i > maxIter %| ((priorDelta - delta)/priorDelta) < tolerance.delta 
       return
   end
   
   regridBool = computeJacobian(gridObject, U, delta, pertubation, tolerance);
   if regridBool == "false"
        U.x = U.x + pertubation.x .* delta;
        U.y = U.y + pertubation.y .* delta;
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
  
    visualize(U.x, U.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate);
    disp("Displaying Displacement Vector fields");
    pause(2);

    
    % apply the deformation to template
    %[movingRegistered]= imwarp(Template,diss);
    % obtain maximal difference 
    %diffImage = (movingRegistered - gridObject.sampleSource);
    % need to fix this ( find the derivative of the difference with respect
    % to the tolerance
%     [gx, gy] = gradient(double(diffImage));
%     [gxx, gxy] = gradient(gx);
%     [gxy, gyy] = gradient(gy);
    
    pause(2);
    priorDelta = delta;
    i = i +1;
    max(max(U.x))
    max(max(U.y))
end
%%
TemplateOut = 0;
% U.x = U.x ./ tolerance.distanceTolerance;
% U.y = U.y ./ tolerance.distanceTolerance;

for i = 1: length(U.x)
    for j = 1: length(U.y)
        TemplateOut(ceil(abs(x(i) + U.x(i, j))), ceil(abs(y(j) + U.y(i, j)))) = Template(x(i),y(j)); 
    end
end

figure; imagesc(TemplateOut);
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
        SourceGrid(x(i),y(j)) = Source(x(i),y(j));
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
