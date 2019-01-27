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
dir = pwd + "/Data/POPI-DataSet/10/*.dcm";

% perform fluid registration
dataStore = fluidRegistration(dir);


% save to .mat filere
save("registrationResults.mat");
%% Analyze Data
mesh( floor(dataStore(1).x - dataStore(1).displacement.x), floor(dataStore.y - dataStore(1).displacement.y), dataStore(1).x .* 0); grid on; colormap gray;
figure; imshowpair(dataStore(1).deformedImage, dataStore(1).Source, "diff"); colormap gray;

figure; imagesc(dataStore(1).deformedImage);

% intensity Difference between Template and Source
visualizeIntensityDiff3D(dataStore(1).Source - dataStore(1).Template);

% intensity Difference between Source and Deformed Template
visualizeIntensityDiff3D(dataStore(1).Source - dataStore(1).deformedImage);

meshGeneration(dataStore(1).gridData, dataStore(1).displacement, "Template");

%%
x = dataStore(1).x;
y = dataStore(1).x;

testSource = zeros(512, 512);
for i = 1: length(x)
    for j = 1: length(y)
        testSource(x(1,i), y(1,j)) = dataStore(1).Source(x(1,i), y(1,j));    
    end
end
visualizeIntensityDiff3D(dataStore(1).deformedImage - testSource);

