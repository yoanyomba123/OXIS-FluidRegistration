function mit18086_navierstokes
%MIT18086_NAVIERSTOKES
%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

% 07/2007 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.
%-----------------------------------------------------------------------
Re = 1e2;     % Reynolds number
dt = 1e0;    % time step
tf = 4e2;    % final time
lx = 512;       % width of box
ly = 512;       % height of box
nx = 100;      % number of x-gridpoints
ny = 100;      % number of y-gridpoints
nsteps = 10;  % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx-1,ny); V = zeros(nx,ny-1);
% boundary conditions
uN = x*0+1;    vN = avg(x)*0;
uS = x*0;      vS = avg(x)*0;
uW = avg(y)*0; vW = y*0;
uE = avg(y)*0; vE = y*0;
%-----------------------------------------------------------------------
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
      [uW;zeros(nx-3,ny);uE]/hy^2);
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);

fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); 
Rp = chol(Lp(perp,perp)); 
Rpt = Rp';
Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:nt
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   
   % extend the velocity matrices by rows and columns of boundary points
   Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   
   % compute the partials in the following manner (UV)x and (UV)y
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg(Ue(:,2:end-1));   Ud = diff(Ue(:,2:end-1))/2;
   Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
   
   % compute the partials in the following manner (U^2)x and (U^2)y
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   
   % update interior points 
   U = U-dt*(UVy(2:end-1,:)+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);
   rhs = reshape(V+Vbc,[],1);
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy;
%    
   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
   if k==1|floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
      % stream function
      rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
      q(perq) = Rq\(Rqt\rhs(perq));
      Q = zeros(nx+1,ny+1);
      Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
      clf, contourf(avg(x),avg(y),Ua,20,'w-'), hold on
      contour(x,y,Q',20,'k-');
      Ue = [uS' avg([uW;U;uE]')' uN'];
      Ve = [vW;avg([vS' V vN']);vE];
      Len = sqrt(Ue.^2+Ve.^2+eps);
      quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-')
      hold off, axis equal, axis([0 lx 0 ly])
      %p = sort(p); caxis(p([8 end-7]))
      title(sprintf('Re = %0.1g   t = %0.2g',Re,k*dt))
      drawnow
   end
end
fprintf('\n')

%=======================================================================

function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;





%     % compute the gradient of the divergence
%     % Centred finite difference matrices
%     delta_x_N = 1/(2*Dx)*(diag(ones(N-1,1),1) - diag(ones(N-1,1),-1));
%     delta_y_M = 1/(2*Dy)*(diag(ones(M-1,1),1) - diag(ones(M-1,1),-1));
% 
%     % Cast U as a vector
%     U_vec = U_mat(:);
% 
%     % Mixed derivative operator
%     A = kron(delta_y_M,delta_x_N);
% 
%     U_xy_num = A*U_vec;
%     U_xy_matrix = reshape(U_xy_num,N,M);
% 
%     % treat nonlinear terms
%     gamma = min(1.2*dt*max(max(max(abs(Vx)))/dx,max(max(abs(Vy)))/dy),1);
%     
%     % extend the velocity matrices by rows and columns of boundary points
%     Vxe = [uW;Vx;uE]; Vxe = [2*uS'-Vxe(:,1) Vxe 2*uN'-Vxe(:,end)];
%     Vye = [vS' Vy vN']; Vye = [2*vW-Vye(1,:);Vye;2*vE-Vye(end,:)];
%     % Treat Nonlinear Terms
%     constant = (mu + lambda);
%     
%     nlVx = avg(Vxe')'; nlVdx = diff(Vxe')';
%     nlVy  = avg(Vye); nlVdy = diff(Vye);
%     
%     
%     
%     %nonLinVx = diff((nlVx .* nlVdx) + (nlVy .* nlVdx))/dx;
%     %nonlineVy = diff((nlVx .* nlVdy) + (nlVy .* nlVdy))/dy;
%     
%     nonLinVx = diff(nlVx.*nlVy-gamma*abs(nlVx).*nlVdy)/dx;
%     nonLinVy = diff((nlVx.*nlVy-gamma*nlVdx.*abs(nlVy))')'/dy;
%    
%     nlVx = avg(Vxe(:,2:end-1));   nlVdx = diff(Vxe(:,2:end-1));
%     nlVy = avg(Vye(2:end-1,:)')'; nlVdy = diff(Vye(2:end-1,:)')';
%     
%     % compute the partials in the following manner (U^2)x and (U^2)y
%     nonLinV2x = diff(nlVx.^2-gamma*abs(nlVx).*nlVdx)/dx;
%     nonLinV2y = diff((nlVy.^2-gamma*abs(nlVy).*nlVdy)')'/dy;
%     
%     % update interior points 
%     Vx = Vx-dt*(nonLinVy(2:end-1,:)+nonLinV2x);
%     Vy = Vy-dt*(nonLinVx(:,2:end-1)+nonLinV2y);
%     
% %      % implicit viscosity
% %      rhs = reshape(Vx+Vxbc,[],1);
% %      u(peru) = Ru\(Rut\rhs(peru));
% %      Vx = reshape(u,numPointsX-1,numPointsY);
% %      rhs = reshape(Vy+Vybc,[],1);
% %      v(perv) = Rv\(Rvt\rhs(perv));
% %      Vy = reshape(v,numPointsX,numPointsY-1);
%%
% clc; clear all; close all;
% clc 
% clear
% myu=0.000012;
% den=1.29;
% x0=0;
% y0=0;
% LX=1;
% LY=1;
% M=20;
% N=20;
% dx=(x0+LX)/M;
% dy=(y0+LY)/N;
% [x,y]=meshgrid(x0:dx:LX,y0:dy:LY);
% plot(x,y,'*r');hold on;grid on
% [xx,yy]=meshgrid(0.1:0.1:1.1,0.1:0.1:1.1);
% plot(xx,yy,'*k');
% 
% for i=1:M;
%     for j=1:N;
%         u(i,j)=randn(1,1);
%         v(i,j)=randn(1,1);
%         uu(i,j)=3;
%         vv(i,j)=1;
%         U(i,j)=u(i,j)+uu(i,j);
%         V(i,j)=v(i,j)+vv(i,j);
%     end
% end
% 
% for i=1:M;
%     for j=1:N;
%         x1(i,j)=i*dx;
%         y1(i,j)=j*dy;
%     end
% end
% 
% for i=1:M;
%     for j=1:N;
%         dudy(i,j)=(u(i+1)-u(i))/dy;
%         dvdx(i,j)=(v(i+1)-v(i))/dx;
%         duudy(i,j)=(uu(i+1)-uu(i))/dy;
%         dvvdx(i,j)=(vv(i+1)-vv(i))/dx;
%         s(i,j)=0.5*(dudy(i,j)+dvdx(i,j));
%         ss(i,j)=0.5*(duudy(i,j)+dvvdx(i,j));
%     end
% end
% 
% for i=1:M;
%     for j=1:N;
%         if (i==j)
%             c=1;
%         elseif (i~=j)
%             c=0;
%         end
%         px(i,j)=0.1*c*randn(1,1);
%         py(i,j)=0.1*c*randn(1,1);
%         ppx(i,j)=0.3*c;
%         ppy(i,j)=0.1*c;
%         Px(i,j)=px(i,j)+ppx(i,j);
%         Py(i,j)=py(i,j)+ppy(i,j);
%         P(i,j)=Px(i,j)+Py(i,j);
%     end
% end
% 
% for i=1:M;
%     for j=1:N;
%         Mstress(i,j)=-(ppx(i,j)/dx+ppy(i,j)/dy)+2*myu*ss(i,j);
%         stressf(i,j)=-(px(i,j)/dx+py(i,j)/dy)+2*myu*s(i,j);
%         for a=1:10;
%         uuu(a)=randn(1,1);
%         vvv(a)=randn(1,1);
%         end
%         um(i,j)=mean(uuu);
%         vm(i,j)=mean(vvv);
%         rt=um(i,j)*vm(i,j);
%         rtt=mean(rt);
%         rrtt=den*rtt;
%         rrrtt=mean(rrtt);
%         ReynoldsStress(i,j)=den*rrrtt;
%     end
% end
% 
% [C,h] = contour(x1,y1,P);
% set(1, 'units', 'centimeters', 'pos', [0 0 120.5 100])
% colormap summer
% 
% quiver(x1,y1,U,V)
