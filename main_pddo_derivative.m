% 增加路径
addpath('PD_Model_Generate_Version2.0');
addpath('PDDO_LIB');
sL =[0.14,0.14]*5;
ndiv =  [7,7];
[coor, pv, ih, jh, ap] = regFastNewHorizon(sL,ndiv, 3.015, 0);
dx = max(sL./ndiv);
mvalue = 3.015;
delta = mvalue*dx;
bl0 = vecnorm(coor(jh,:)-coor(ih,:),2,2);
wf = exp(-bl0/delta); %权函数
vc = (bl0<=delta-dx/2).*1 + (bl0>delta-dx/2).*(delta - bl0 + dx/2)/dx;
wfvj = wf.*vc.*pv(jh); % 权函数乘以体积
clear wf vc
ndim  = size(coor,2);
nSeries = 2; 
gf = getPDDOFun(ndim, nSeries, coor, ih, jh ,ap, wfvj);
g10 = gf(:,1).*wfvj;
g01 = gf(:,2).*wfvj;
g20 = gf(:,3).*wfvj;
g11 = gf(:,4).*wfvj;
g02 = gf(:,5).*wfvj;

epsilon = 1e-2;
f = @(x,y) epsilon*x + epsilon*y + epsilon*x.*y + epsilon*x.^2 + epsilon*x.^3;
dfdx = @(x,y) epsilon + epsilon*y + 2*epsilon*x + 3*epsilon*x.^2;
dfdx2 = @(x,y) 2*epsilon*ones(size(x)) + 6*epsilon*x;
x = coor(:,1);
y = coor(:,2);
fv = f(x,y);
dfv = fv(jh) - fv(ih);
dfdx_cal  = accumarray(ih, dfv.*g10);
dfdx2_cal = accumarray(ih, dfv.*g20);
dfdx_real = dfdx(x,y);
dfdx2_real = dfdx2(x,y);
x = reshape(x,ndiv);
y = reshape(y,ndiv);
dfdx_cal = reshape(dfdx_cal,ndiv);
dfdx2_cal = reshape(dfdx2_cal,ndiv);
dfdx_real = reshape(dfdx_real,ndiv);
dfdx2_real = reshape(dfdx2_real,ndiv);
figure(1); clf
subplot(2,3,1)
pcolor(x,y,dfdx_cal); colormap(jet); colorbar('southoutside'); axis equal; title('Cal df/dx')
axis off
subplot(2,3,2)
pcolor(x,y,dfdx_real); colormap(jet); colorbar('southoutside'); axis equal; title('Real df/dx')
axis off
subplot(2,3,3)
pcolor(x,y,dfdx_cal - dfdx_real); colormap(jet); colorbar('southoutside'); axis equal; title('Error')
axis off

subplot(2,3,4)
pcolor(x,y,dfdx2_cal); colormap(jet); colorbar('southoutside'); axis equal; title('Cal d^2f/dx^2')
axis off
subplot(2,3,5)
pcolor(x,y,dfdx2_real); colormap(jet); colorbar('southoutside'); axis equal; title('Real d^2f/dx^2')
axis off
subplot(2,3,6)
pcolor(x,y,dfdx2_cal - dfdx2_real); colormap(jet); colorbar('southoutside'); axis equal; title('Error')
axis off