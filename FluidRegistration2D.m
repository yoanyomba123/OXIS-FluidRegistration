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
templateSourceDiff = Template - Source;
 
% display and visualize both images
figure; imagesc([Template, Source]);

% display and visualize difference between both images
figure; imshowpair(Template, Source, "diff");

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
tolerance.deformationTolerance = 0.5;
tolerance.jacobainTolerance = 0.025;
tolerance.distanceTolerance = 1e-5;

maxIter = 20;

% grid definition
[rows, cols] = size(Template);
gridObject = struct();
gridObject.numXPoints = ceil(rows/10);
gridObject.numYPoints = ceil(cols/10);
gridObject.grid = struct();

% generate points that are not on the boundary of the image
x = linspace(0, rows-2, gridObject.numXPoints); x = ceil(x);
y = linspace(0, cols-2, gridObject.numYPoints); x = ceil(y);
[X,Y] = meshgrid(x, y);
gridObject.rows = rows;
gridObject.cols = cols;
gridObject.x = x;
gridObject.y = y;
gridObject.grid.x = X;
gridObject.grid.y = Y;
gridObject.dx = ceil(gridObject.rows/gridObject.numXPoints);
gridObject.dy = ceil(gridObject.cols/gridObject.numYPoints);

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
%% Obtain True Template
sampleTemplate = gridObject.grid.x;
for i = 2:length(x)
    for j = 2:length(y)
        sampleTemplate(i,j) = Template(x(i), y(j));
        sampleSource(i,j) = Source(x(i), y(j));
    end
end
gridObject.sampleTemplate = sampleTemplate;
gridObject.sampleSource = sampleSource;
%% Algorithm
wk = struct();
for i = 1:1;
    % perform 2D interpolation
    % turn this into a function
    wRegrid.x = interp2(yQ{regridCounter}.x, gridObject.grid.x - U.x);
    wRegrid.y = interp2(yQ{regridCounter}.y, gridObject.grid.y - U.y);
    wK{i} = wRegrid;
    
    % turn this into a function
    tRegrid.x = interp2(gridObject.sampleTemplate, gridObject.grid.x - wRegrid.x - U.x);
    tRegrid.y = interp2(gridObject.sampleTemplate, gridObject.grid.y - wRegrid.y - U.y);
    tK{i} = tRegrid;
    
    if(i > 1)
       if((gridObject.sampleSource - tk{i-1}.x) <=  (gridObject.sampleSource - tk{i}.x).*tolerance.distanceTolerance && ...
               (gridObject.sampleSource - tk{i-1}.x) <=  (gridObject.sampleSource - tk{i}.x).*tolerance.distanceTolerance)
            exit;
       end
    end
    
    force = forceField(Template, Source, U, gridObject, "none");
    
    drawnow
    % visualize the force field on the image
    visualize(force.x, force.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate - gridObject.sampleSource);
    disp("Displaying force fields");
    pause(3);
    
    % obtain V
    % Turn this into a function
    Sx = stencil.S11 + stencil.S12;
    Sy = stencil.S21 + stencil.S22;
    
    A.x = conv2(gridObject.sampleTemplate, Sx, 'same');
    A.y = conv2(gridObject.sampleTemplate, Sy, 'same');
    A.FFTx = fftMatOperator .* A.x .* fftMatInvOperator;
    A.FFTy = fftMatOperator .* A.y .* fftMatInvOperator;
    
    Dx = pinv(A.FFTx);
    Dy = pinv(A.FFTy);
    
    V.x = Dx .* force.x;
    V.y = Dy .* force.y;
    % visualize the velocity field on the image
    visualize(V.x, V.y, gridObject.grid.x, gridObject.grid.y, gridObject.sampleTemplate - gridObject.sampleSource);
    disp("Displaying Velocity Vector fields");
    pause(3);
    
end