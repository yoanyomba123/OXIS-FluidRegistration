%% Viscous Fluid Image Registrion
% Author: D Yoan L. Mekontchou Yomba
% Date: 11/16/2018
% Purpose:  
%   The purpose of this script is to implement a fluid registration model
%   in matlab. The Registration model makes use of  the Eularian reference
%   frame and uses the sum of square difference as a cost function.
%% Clear Up Workspace
clc; clear; close all;

%% Add all paths to current workspace recursively
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Perform File Processing
dir = pwd + "/Data/POPI-DataSet/00/*.dcm";
files = readImagesFromDirectory(dir);

numpoints = 100;
iter = 500;
mu = 10; lambda = 10;

% tolerance definition
tolerance = struct();
tolerance.deformationTolerance = 50;
tolerance.jacobianTolerance = 0.5;
tolerance.distanceTolerance = 1e-2;
tolerance.mse = 1e-12;

% read all images and augment them then register
for i = 1:2
    % transform / warp image to create template
    Template = imrotate(files{i},-60,'bilinear','crop'); % rotate the template image

    % obtain source image
    Source = files{i};
    
    % perform registration of images
    [Template, Source, preMSE, postMSE, preIntDiff, postIntDiff, TemplateOut] = fluidRegistration2DImpl(Template, Source, numpoints, iter, mu, lambda, tolerance);
    
    % store relevant stats in a struct for each image
    % mse pre and post registration, intensity difference pre and post
    % registration, deformed image, as well as template and source images
    regRes = struct();
    regRes(i).Template = Template;
    regRes(i).Source = Source;
    regRes(i).preRegMSE = preMSE;
    regRes(i).postRegMSE = postMSE;
    regRes(i).preIntDiff = preIntDiff;
    regRes(i).postIntDiff = postIntDiff;
    regRes(i).deformedImage = TemplateOut;
end

%% Visualize Deformed Image
figure; imagesc(regRes(2).deformedImage);
