function [coor, pv, ih, jh, accumPos] = regFastNewHorizon(sideLength, ndiv, mvalue, isboundary)
%% 一种比较快速的计算规则排布物质点的邻域的算法
%% 作者：万冀
%% 单位：武汉大学
%% 版本更新信息
%% 2020/11/06 Version 2.0 更新了邻域元胞转换成邻域向量的算法
%% 2021/05/05 Version 2.0 纠正了邻域物质点的指示方法
ndim = numel(sideLength);
dx = sideLength./ndiv;
delta = mvalue * max(dx);
if(ndim==3) % 三维坐标一般先是z从0-nz再是y从0-ny最后x从0-nx 
    if(isboundary==0)
        x = (1-ndiv(1):2:ndiv(1)-1)*0.5*dx(1);
        y = (1-ndiv(2):2:ndiv(2)-1)*0.5*dx(2);
        z = (1-ndiv(3):2:ndiv(3)-1)*0.5*dx(3); 
    elseif(isboundary~=0)
        x = (-ndiv(1):2:ndiv(1))*0.5*dx(1);
        y = (-ndiv(2):2:ndiv(2))*0.5*dx(2);
        z = (-ndiv(3):2:ndiv(3))*0.5*dx(3);
    end
    if(isboundary~=0)
        ndiv = ndiv+1;
    end
    np = prod(ndiv);
    [Y,X,Z] = meshgrid(y,x,z);
    coor = [reshape(X,np,1), reshape(Y,np,1),reshape(Z,np,1)];
    [J, I, K] = meshgrid(1:ndiv(2),1:ndiv(1), 1:ndiv(3));
    N = sub2ind(ndiv, I, J, K);
    endiv = floor(delta./dx);
    p = endiv(1);
    q = endiv(2);
    r = endiv(3);
    [j,i,k] = meshgrid(-q:q, -p:p, -r:r);
    idx = (i.^2*dx(1)^2+j.^2*dx(2)^2+k.^2*dx(3)^2)<=delta^2 ;
    idx(p+1,q+1,r+1) = false;
    horizon = cell(ndiv);
    for i = 1:1:ndiv(1)
        for j= 1:1:ndiv(2)
            for k = 1:1:ndiv(3)
                p0 = max(i-p,1):min(i+p, ndiv(1));
                q0 = max(j-q,1):min(j+q, ndiv(2));
                r0 = max(k-r,1):min(k+r, ndiv(3));
                eN = N(p0,q0,r0);
                padd = idx(p0-i+p+1,q0-j+q+1,r0-k+r+1);
                horizon{i,j,k} = eN(padd);
            end
        end
    end
elseif(ndim==2) % y从0-ny最后x从0-nx
    if(isboundary==0)
        x = (1-ndiv(1):2:ndiv(1)-1)*0.5*dx(1);
        y = (1-ndiv(2):2:ndiv(2)-1)*0.5*dx(2);
    elseif(isboundary~=0)
        x = (-ndiv(1):2:ndiv(1))*0.5*dx(1);
        y = (-ndiv(2):2:ndiv(2))*0.5*dx(2);
    end
    if(isboundary~=0)
        ndiv = ndiv+1;
    end
    np = prod(ndiv);
    [Y,X] = meshgrid(y,x);
    coor = [reshape(X,np,1), reshape(Y,np,1)];
    [J, I] = meshgrid(1:ndiv(2),1:ndiv(1));
    N = sub2ind(ndiv, I, J);
    endiv = floor(delta./dx);
    p = endiv(1);
    q = endiv(2);
    [j,i] = meshgrid(-q:q, -p:p);
    idx = (i.^2*dx(1)^2+j.^2*dx(2)^2)<=delta^2 ;
    idx(p+1,q+1) = false;
    horizon = cell(ndiv);
    for i = 1:1:ndiv(1)
        for j= 1:1:ndiv(2)
                p0 = max(i-p,1):min(i+p, ndiv(1));
                q0 = max(j-q,1):min(j+q, ndiv(2));
                eN = N(p0,q0);
                padd = idx(p0-i+p+1,q0-j+q+1);
                horizon{i,j} = eN(padd);
        end
    end
elseif(ndim==1)
    if(isboundary==0)
        x = (1-ndiv(1):2:ndiv(1)-1)*0.5*dx(1);
    elseif(isboundary~=0)
        x = (-ndiv(1):2:ndiv(1))*0.5*dx(1);
    end
    coor = x';
    if(isboundary~=0)
        ndiv = ndiv+1;
    end
    np = ndiv;
    endiv = floor(delta./dx);
    p = endiv(1);
    horizon = cell(ndiv(1),1);
    for i = 1:1:ndiv(1)
                p0 = [max(i-p,1):i-1,i+1:min(i+p, ndiv(1))]';
                horizon{i} = p0;
    end
end
pv = ones(np,1)*prod(dx);
horizon = reshape(horizon, np,1);
npph = zeros(np,1);
for i = 1:1:np
    npph(i) = numel(horizon{i});
end
accumPos = [0; cumsum(npph)];
ih = zeros(accumPos(end),1); jh=ih;
for i = 1:1:np
    jh(accumPos(i)+1:accumPos(i+1)) = horizon{i,1};
    ih(accumPos(i)+1:accumPos(i+1)) = i;
end
end