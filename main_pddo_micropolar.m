clc;clear
% 增加路径
addpath('PD_Model_Generate_Version2.0');
addpath('PDDO_LIB');

fileDir = 'resFolder2\';

%% displacement load
applyU1U2arr = 1e-4*[0:2:20, 20+0.02:0.02:50];
crackLen = 0.5d0;
nt = numel(applyU1U2arr);

%% material constants
lb = 1*0.5;
couplingNumber = 0.50;
scr = 0.01;
pr = 0.2;
emod = 12e6;
smod = emod/(2*(1+pr));
gamma = 4*smod*lb^2;
kappa = 2*couplingNumber^2/(1-couplingNumber^2)*smod;
mu = smod - kappa/2;
lambda = pr*emod/(1+pr)/(1-2*pr);
alpha = 0;
pd_model.kappa = kappa;
pd_model.mu = mu;
pd_model.lb = lb;
pd_model.lambda = lambda;
pd_model.nt = nt;
pd_model.crackLen = crackLen;
pd_model.coupling_number = couplingNumber;
pd_model.pr = pr;
pd_model.emod = emod;
pd_model.scr = scr;

if(~exist(fileDir, 'dir'))
    mkdir(fileDir)
end
sL = [1,1]; % side length
ndiv = [101, 101];
mvalue = 3.015;
isboundaryIncluded = 1;
[coor, pv, ih, jh, ap] = regFastNewHorizon(sL, ndiv, mvalue, isboundaryIncluded);

pd_model.coor = coor;
pd_model.pv  = pv;
pd_model.ih  = ih;
pd_model.jh = jh;
start_point = [-sL(1)*0.5-eps,0.0+eps];
end_point = [-sL(1)*0.5+crackLen+eps,0.0+eps];
fail = preCrack2D(coor, ih, jh, start_point, end_point); % prefabricated crack
dmg = accumarray(ih, pv(jh).*fail)./accumarray(ih, pv(jh));
scatter(coor(:,1), coor(:,2), 20, dmg, 'filled')
axis equal
axis off
colormap(jet)
colorbar
dx = max(sL./ndiv);
delta = mvalue*dx;
bl0 = vecnorm(coor(jh,:)-coor(ih,:),2,2);
wf = exp(-bl0/delta); %权函数
vc = (bl0<=delta-dx/2).*1 + (bl0>delta-dx/2).*(delta - bl0 + dx/2)/dx;
wfvj = wf.*vc.*pv(jh); % 权函数乘以体积
clear wf vc
ndim  = size(coor,2);
nSeries = 2;
% For [ndim, nSeries] = [2,2] Case, [10][01][20][11][02] 是对的
gf = getPDDOFun(ndim, nSeries, coor, ih, jh ,ap, wfvj.*~fail);

downNode = find(coor(:,2)<-sL(2)/2+eps);
upNode = find(coor(:,2)>sL(2)/2-eps);
leftNode = find(coor(:,1)<-sL(1)/2+eps&(~(coor(:,2)<-sL(2)/2+eps))&(~(coor(:,2)>sL(2)/2-eps)));
leftNodeUp = leftNode(coor(leftNode,2)>0+eps);
leftNodeDown = leftNode(coor(leftNode,2)<=0+eps);
rightNode = find(coor(:,1)>sL(1)/2-eps&(~(coor(:,2)<-sL(2)/2+eps))&(~(coor(:,2)>sL(2)/2-eps)));
crackNode = find(coor(:,1)> -sL(1)/2+eps & coor(:,1)< -sL(1)/2+crackLen-eps&...
    coor(:,2)<dx+eps&coor(:,2)>-dx-eps);
crackNodeUp = crackNode(coor(crackNode,2)>0+eps);
crackNodeDown = crackNode(coor(crackNode,2)<=0+eps);
scatter(coor(downNode,1),coor(downNode,2),10,'r','filled');
hold on
scatter(coor(upNode,1),coor(upNode,2),10,'b','filled');
scatter(coor(leftNodeUp,1),coor(leftNodeUp,2),10,'m','filled');
scatter(coor(leftNodeDown,1),coor(leftNodeDown,2),10,'k','filled');
scatter(coor(rightNode,1),coor(rightNode,2),10,'c','filled');
%scatter(coor(crackNode,1),coor(crackNode,2),10,'k','filled');
scatter(coor(crackNodeUp,1),coor(crackNodeUp,2),10,'b','filled');
scatter(coor(crackNodeDown,1),coor(crackNodeDown,2),10,'r','filled');


nDownNode = numel(downNode);
nUpNode = numel(upNode);
nLeftNode = numel(leftNode);
nRightNode = numel(rightNode);

np = size(coor,1); % material point number

X = coor(jh,:)-coor(ih,:); % relative position

save ([fileDir,'pd_model'], 'pd_model')
for tt = 1:1:nt
    %% [ndim, nSeries]=[2,2] => [10][01][20][11][02]
    g10 = gf(:,1).*wfvj.*~fail;
    g01 = gf(:,2).*wfvj.*~fail;
    g20 = gf(:,3).*wfvj.*~fail;
    g11 = gf(:,4).*wfvj.*~fail;
    g02 = gf(:,5).*wfvj.*~fail;
    
    
    
    npa = (1:np)';
    p11 = (lambda+mu)*g20 + (mu+kappa)*(g20+g02);
    p12 = (lambda+mu)*g11;
    p13 = kappa*g01;
    p21 = (lambda+mu)*g11;
    p22 = (lambda+mu)*g02 + (mu+kappa)*(g20+g02);
    p23 = -kappa*g10;
    p31 = -kappa*g01;
    p32 = +kappa*g10;
    p33 = gamma*(g20+g02);
    Kt = [
        3*ih-2, 3*jh-2, p11;
        3*npa-2, 3*npa-2, -accumarray(ih, p11);
        3*ih-2, 3*jh-1, p12;
        3*npa-2, 3*npa-1, -accumarray(ih, p12);
        3*ih-2, 3*jh,   p13;
        3*npa-2, 3*npa,   -accumarray(ih, p13);
        3*ih-1, 3*jh-2, p21;
        3*npa-1, 3*npa-2, -accumarray(ih, p21);
        3*ih-1, 3*jh-1, p22;
        3*npa-1, 3*npa-1, -accumarray(ih, p22);
        3*ih-1, 3*jh,   p23;
        3*npa-1, 3*npa,   -accumarray(ih, p23);
        3*ih, 3*jh-2,   p31;
        3*npa, 3*npa-2,   -accumarray(ih, p31);
        3*ih, 3*jh-1,   p32;
        3*npa, 3*npa-1,   -accumarray(ih, p32);
        3*ih, 3*jh,     p33;
        3*npa, 3*npa,     -accumarray(ih, p33)-2*kappa;
        ];
    totDOF = np*3;
    Ktot = sparse(Kt(:,1),Kt(:,2),Kt(:,3),totDOF, totDOF);
    % (Ktot)
    
    %% 各个边界物质点数目
    nleftNodeUp = numel(leftNodeUp);
    nleftNodeDown = numel(leftNodeDown);
    ncrackNodeUp = numel(crackNodeUp );
    ncrackNodeDown = numel(crackNodeDown);
    nupNode = numel(upNode);
    ndownNode = numel(downNode);
    nrightNode = numel(rightNode);
    
    nBound = 3*(nleftNodeUp+nleftNodeDown+ncrackNodeUp+ncrackNodeDown+nupNode+ndownNode+nrightNode);
    nPos = 0;
    %% 初始化边界条件矩阵和力矩阵
    G = sparse([],[],[], nBound, totDOF);
    F = sparse([],[],[],nBound,1);
    %% 施加位移荷载
    apply_disp_u1 = 0;
    apply_disp_u2 = applyU1U2arr(tt);
    F = F + sparse((nleftNodeUp+nleftNodeDown+ncrackNodeUp+ncrackNodeDown+ndownNode)*3+2:3:3*...
        (nleftNodeUp+nleftNodeDown+ncrackNodeUp+ncrackNodeDown+ndownNode+nupNode),1,apply_disp_u2,nBound,1);
    nCount = 0;
    for k = 1:1:numel(leftNodeUp)
        i = leftNodeUp(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = -1; n2 = 0;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap); n1*lambda*g01(eap)+n2*mu*g10(eap); ...
            -sum(n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap)); ...
            -sum(n1*lambda*g01(eap)+n2*mu*g10(eap));...
            -n1*kappa],...
            nBound, totDOF)+... % t1 = 0;
            sparse( nCount*3-1, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*mu*g01(eap)+n2*lambda*g10(eap); n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap); ...
            -sum(n1*mu*g01(eap)+n2*lambda*g10(eap));...
            -sum(n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap));...
            n2*kappa],...
            nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    for k = 1:1:numel(leftNodeDown)
        i = leftNodeDown(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = -1; n2 = 0;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap); n1*lambda*g01(eap)+n2*mu*g10(eap); ...
            -sum(n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap)); ...
            -sum(n1*lambda*g01(eap)+n2*mu*g10(eap));...
            -n1*kappa],...
            nBound, totDOF)+... % t1 = 0;
            sparse( nCount*3-1, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*mu*g01(eap)+n2*lambda*g10(eap); n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap); ...
            -sum(n1*mu*g01(eap)+n2*lambda*g10(eap));...
            -sum(n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap));...
            n2*kappa],...
            nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    
    for k = 1:1:numel(crackNodeUp)
        i = crackNodeUp(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = 0; n2 = -1;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap); n1*lambda*g01(eap)+n2*mu*g10(eap); ...
            -sum(n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap)); ...
            -sum(n1*lambda*g01(eap)+n2*mu*g10(eap));...
            -n1*kappa],...
            nBound, totDOF)+... % t1 = 0;
            sparse( nCount*3-1, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*mu*g01(eap)+n2*lambda*g10(eap); n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap); ...
            -sum(n1*mu*g01(eap)+n2*lambda*g10(eap));...
            -sum(n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap));...
            n2*kappa],...
            nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    for k = 1:1:numel(crackNodeDown)
        i = crackNodeDown(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = 0; n2 = 1;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap); n1*lambda*g01(eap)+n2*mu*g10(eap); ...
            -sum(n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap)); ...
            -sum(n1*lambda*g01(eap)+n2*mu*g10(eap));...
            -n1*kappa],...
            nBound, totDOF)+... % t1 = 0;
            sparse( nCount*3-1, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*mu*g01(eap)+n2*lambda*g10(eap); n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap); ...
            -sum(n1*mu*g01(eap)+n2*lambda*g10(eap));...
            -sum(n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap));...
            n2*kappa],...
            nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    
    for k = 1:1:numel(downNode)
        i = downNode(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = 0; n2 = -1;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, 3*i-2,1, nBound, totDOF)+...
            sparse( nCount*3-1, 3*i-1, 1, nBound, totDOF)+...
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    for k = 1:1:numel(upNode)
        i = upNode(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = 0; n2 = 1;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, 3*i-2,1, nBound, totDOF)+...
            sparse( nCount*3-1, 3*i-1, 1, nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    
    for k = 1:1:numel(rightNode)
        i = rightNode(k);
        eap = ap(i)+1:ap(i+1);
        eih = ih(eap);
        ejh = jh(eap);
        n1 = 1; n2 = 0;
        nCount = nCount + 1;
        G = G + sparse( nCount*3-2, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap); n1*lambda*g01(eap)+n2*mu*g10(eap); ...
            -sum(n1*(lambda+2*mu+kappa)*g10(eap)+n2*(mu+kappa)*g01(eap)); ...
            -sum(n1*lambda*g01(eap)+n2*mu*g10(eap));...
            -n1*kappa],...
            nBound, totDOF)+... % t1 = 0;
            sparse( nCount*3-1, [ejh*3-2; ejh*3-1; 3*i-2; 3*i-1; 3*i], ...
            [n1*mu*g01(eap)+n2*lambda*g10(eap); n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap); ...
            -sum(n1*mu*g01(eap)+n2*lambda*g10(eap));...
            -sum(n2*(lambda+2*mu+kappa)*g01(eap)+n1*(mu+kappa)*g10(eap));...
            n2*kappa],...
            nBound, totDOF)+... % t2 = 0;
            sparse( nCount*3, 3*i, 1,  nBound, totDOF);
    end
    
    U = [Ktot, G'; G, sparse([],[],[],nBound,nBound)]\[sparse([],[],[],totDOF,1);F];
    u1 = U(1:3:totDOF);
    u2 = U(2:3:totDOF);
    omega3 = U(3:3:totDOF);
    
    s = (vecnorm(X+[u1(jh)-u1(ih)+X(:,2).*(omega3(ih)+omega3(jh))*0.5, ...
        u2(jh)-u2(ih)-X(:,1).*(omega3(ih)+omega3(jh))*0.5],2,2) - bl0)./bl0;
    %s = (vecnorm([X, u3(jh)-u3(ih)],2,2) - bl0)./bl0;
    fprintf('tt=%d\n',tt)
    fprintf('最大剩余键键长：%f\n',max(s(~fail)));
    failbefore = fail;
    fail = fail | (s>scr);
    nbrokenBond = sum(fail~=failbefore);
    fprintf('新增断键个数%d\n',nbrokenBond)
    dmg = accumarray(ih, fail.*wfvj)./accumarray(ih, wfvj) ;
    s11 = accumarray(ih, (lambda+2*mu+kappa)*g10.*(u1(jh)-u1(ih)) + lambda*g01.*(u2(jh)-u2(ih)));
    s12 = accumarray(ih, (mu+kappa)*g10.*(u2(jh)-u2(ih))+mu*g01.*(u1(jh)-u1(ih))) - kappa*omega3;
    s21 = accumarray(ih, (mu+kappa)*g01.*(u1(jh)-u1(ih))+mu*g10.*(u2(jh)-u2(ih))) + kappa*omega3;
    s22 = accumarray(ih, (lambda+2*mu+kappa)*g01.*(u2(jh)-u2(ih)) + lambda*g10.*(u1(jh)-u1(ih)));
    m31 = accumarray(ih, gamma*(omega3(jh)-omega3(ih)).*g10);
    m32 = accumarray(ih, gamma*(omega3(jh)-omega3(ih)).*g01);
    result.u1 = u1;
    result.u2 = u2;
    result.omega3 = omega3;
    result.dmg = dmg;
    result.apply_disp_u1 = apply_disp_u1;
    result.apply_disp_u2 = apply_disp_u2;
    result.s11 = s11;
    result.s12 = s12;
    result.s21 = s21;
    result.s22 = s22;
    result.m31 = m31;
    result.m32 = m32;
    if ( mod(tt,10)==0 ) % output result every 10 step
        save ([fileDir,'result_step_',num2str(tt)], 'result')
    end
    if(mod(tt,10)==0)
        figure(108)
        clf
        scatter(coor(:,1),coor(:,2),5,dmg,'filled')
        axis equal
        title(['dmg:u1=',num2str(apply_disp_u1*1e3),'mm lc=',num2str(lb)])
        colormap(jet)
        colorbar
        pause(0.000001)
    end
end

%% some post processing code
% figure(1)
% clf
% scatter(coor(:,1),coor(:,2),10,omega3,'filled')
% axis equal
% colormap(jet)
% colorbar
%
% figure(2)
% clf
% scatter(coor(:,1),coor(:,2),10,u1,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('u_1')
%
% figure(3)
% clf
% scatter(coor(:,1),coor(:,2),10,u2,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('u_2')
%
% s11 = accumarray(ih, (lambda+2*mu+kappa)*g10.*(u1(jh)-u1(ih)) + lambda*g01.*(u2(jh)-u2(ih)));
% s12 = accumarray(ih, (mu+kappa)*g10.*(u2(jh)-u2(ih))+mu*g01.*(u1(jh)-u1(ih))) - kappa*omega3;
% s21 = accumarray(ih, (mu+kappa)*g01.*(u1(jh)-u1(ih))+mu*g10.*(u2(jh)-u2(ih))) + kappa*omega3;
% s22 = accumarray(ih, (lambda+2*mu+kappa)*g01.*(u2(jh)-u2(ih)) + lambda*g10.*(u1(jh)-u1(ih)));
%
% m31 = accumarray(ih, gamma*(omega3(jh)-omega3(ih)).*g10);
% m32 = accumarray(ih, gamma*(omega3(jh)-omega3(ih)).*g01);
%
%
% figure(30)
% clf
% scatter(coor(:,1),coor(:,2),10,s11,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Stress \sigma_{11}(Pa)')
%
% figure(31)
% clf
% clf
% scatter(coor(:,1),coor(:,2),10,s12,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Stress \sigma_{12}(Pa)')
%
% figure(32)
% clf
% scatter(coor(:,1),coor(:,2),10,s21,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Stress \sigma_{21}(Pa)')
%
% figure(33)
% clf
% scatter(coor(:,1),coor(:,2),10,s22,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Stress \sigma_{22}(Pa)')
%
% figure(34)
% clf
% scatter(coor(:,1),coor(:,2),10,m31,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Couple stress m_{31}(Pa·m)')
%
% figure(35)
% clf
% scatter(coor(:,1),coor(:,2),10,m32,'filled')
% axis equal
% colormap(jet)
% colorbar
% title('Couple stress m_{32}(Pa·m)')