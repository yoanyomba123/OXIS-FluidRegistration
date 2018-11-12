%% Clear Up Workspace
clc; clear; close all;

% add all paths to current workspace recursively
currentFolder = pwd;
addpath(genpath(currentFolder));

% Read In Images
Template = im2double(dicomread("Data/Template.dcm"));
Source = im2double(dicomread("Data/Source.dcm"));
%% Defining First Set Of Initial Conditions
Re = 1; % Reynolds Number
% Reynolds Number is in this case due in part to the fact
% that the aim is to wholely incorporate the diffusive
% mechaninism into the model without performing alterations
x0=0;
y0=0;

dt = 10; % initialy define the maximal time step
Umax =  1e-1; % define the deformation Limit - Need To come Back
% to this
tInitial = 0;  % define the initial time step
tFinal = 20; % define maximal iterations;
%% Defining Second Set Of Initial Conditions
% grid/mesh size shoud match that of the template 
[rows, cols] = size(Template);
gridLengthX = rows; % grid width
gridlengthY = cols; % grid height
% define number of control points in each direction
numPointsX = 40+1;
numPointsY = 40+1;
% define number of time step until we observe current system
% state visualy
numSteps = 20;
%% Defining Second Set Of Experimental Intial Conditions
numTimeSteps = ceil(tFinal/dt);
dt = ceil(tFinal/numTimeSteps);

% Note that we must include a few points past the border in order
% in order to efficiently perform linear interpolation at image boundaries
% Maked sure to not generate a point for the image edge 

% NOTE TODO: Must Fix this and resize image as well as create a bigger
% range of values in order to adequately perform linear interpolation at
% image edges
x = linspace(0, gridLengthX-2, numPointsX+1);
y = linspace(0, gridlengthY-2, numPointsY+1);

% compute the spacing between each control point in x and y direction
dx = ceil(gridLengthX/numPointsX);
dy = ceil(gridlengthY/numPointsY);

% create a meshgrid
[X, Y] = meshgrid(x,y);
x = ceil(x); y = ceil(y);
%[x,y]=meshgrid(x0:dx:LX,y0:dy:LY);
plot(X,Y,'*r');hold on;grid on
figure; imagesc(Template);
%[xx,yy]=meshgrid(0.1:0.1:1.1,0.1:0.1:1.1);
%% Third Initialization Sequence
% initialize displacement field
Ux = zeros(numPointsX, numPointsY);
Uy = zeros(numPointsX, numPointsY);

% initialize velocity field
Vx = zeros(numPointsX, numPointsY);
Vy = zeros(numPointsX, numPointsY);


%% Solve for velocity based off of the PDE
mu = 1;
lambda = 1;

% boundary conditions
uN = x*0+1;    vN = avg(x)*0;
uS = x*0;      vS = avg(x)*0;
uW = avg(y)*0; vW = y*0;
uE = avg(y)*0; vE = y*0;
%-----------------------------------------------------------------------
Vxbc = dt/Re*([2*uS(2:end-1)' zeros(numPointsX-1,numPointsY-2) 2*uN(2:end-1)']/dx^2+...
      [uW;zeros(numPointsX-3,numPointsY);uE]/dy^2);
Vybc = dt/Re*([vS' zeros(numPointsX,numPointsY-3) vN']/dx^2+...
      [2*vW(2:end-1);zeros(numPointsX-2,numPointsY-1);2*vE(2:end-1)]/dy^2);


fprintf('initialization')
Lp = kron(speye(numPointsY),K1(numPointsX,dx,1))+kron(K1(numPointsY,dy,1),speye(numPointsX));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); 
Rp = chol(Lp(perp,perp)); 
Rpt = Rp';
Lu = speye((numPointsX-1)*numPointsY)+dt/Re*(kron(speye(numPointsY),K1(numPointsX-1,dx,2))+...
     kron(K1(numPointsY,dy,3),speye(numPointsX-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(numPointsX*(numPointsY-1))+dt/Re*(kron(speye(numPointsY-1),K1(numPointsX,dx,3))+...
     kron(K1(numPointsY-1,dy,2),speye(numPointsX)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
Lq = kron(speye(numPointsY-1),K1(numPointsX-1,dx,2))+kron(K1(numPointsY-1,dy,2),speye(numPointsX-1));


centralDiffMat = full(gallery('tridiag',numPointsX,-1,2,-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
for i=1:9;
    % Force Field Computation
    force = forceField(Template, Source, x ,y, Ux, Uy,dx, dy);
    
    % Visualize Force Vector Field At each Location for the Constituent
    % image
    figure; quiver(X, Y, force(:,:,1), force(:,:,2));
    disp("Displaying force fields");
    % compute 2nd order DIFFQ of V wrt xx and yy
    d2Vx_xx = centralDiffMat .* Vx;  d2Vx_yy = centralDiffMat .* Vx'; % multilplying by the transpose gives us differentiation in the y direction
    d2Vy_xx = centralDiffMat .* Vy;  d2Vy_yy = centralDiffMat .* Vy;
    
    % compute 2nd order DIFFQ of V wrt xy
    % Cast U as a vector
    Vx_vec = Vx(:);    Vy_vec = Vy(:);
    % Mixed derivative operator
    Ax = kron(d2Vx_xx,d2Vx_yy);     Ay = kron(d2Vy_xx,d2Vy_yy);
    
    Vx_xy_num = Ax*Vx_vec;    Vy_xy_num = Ay*Vy_vec;

    d2Vx_xy = reshape(Vx_xy_num,numPointsX,numPointsY);
    d2Vy_xy = reshape(Vy_xy_num,numPointsX,numPointsY);

    
    % Compute A
    A11 = ((mu + 2*lambda) .* d2Vx_xx) + (mu .* d2Vx_yy);
    A12 = (mu + lambda) .* d2Vy_xy;
    A21 = (mu + lambda) .* d2Vx_xy;
    A22 = (mu) .* d2Vy_xx + (mu + 2*lambda) .* d2Vy_yy;
    
    % Perform pseudoinverse
    D11 = diag(diag(A11));
    D12 = diag(diag(A12));
    D21 = diag(diag(A21));
    D22 = diag(diag(A22));
    
    Vx = force(1,1,1) .* D11 + force(1,1,2) .*D12;
    Vy = force(1,1,1) .* D21 + force(1,1,2) .*D22;
    figure; quiver(X(2:end, 2:end), Y(2:end, 2:end), Vx, Vy);
    disp("Displaying velocity vector fields \n");
    
    % Compute the pertubation
    Px = Vx -  sum(sum(Vx * diff(Ux)'));
    Py = Vy -  sum(sum(Vy * diff(Uy)'));
    
    % obtain delta based on the quotient of the limit of the deformaton and
    % the pertubation
    ind = Umax/(max( max(max(Px)) ,  max(max(Py))));
    delta = ind - 0.05;
    
    % Update U
    Ux = Ux + delta.*Px;
    Uy = Uy + delta.*Py;
    figure; quiver(X(2:end, 2:end), Y(2:end, 2:end), Ux, Uy);
    disp("Displaying displacement vector fields \n");
    
end

figure; quiver(X(1:end-2, 1:end-1),Y(1:end-2,1:end-1),Vx,Vy');
%%
%% 
clear all; close all;
mit18086_navierstokes()
%%
for i=1:length(X);
    for j = 1:length(Y);
        x0 = [i,j];
        % compute the body field
        F = 7;
        
    end
end







