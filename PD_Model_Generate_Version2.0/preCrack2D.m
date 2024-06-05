function fail = preCrack2D(coor, ih, jh, sp, ep)
%% 预制裂纹
% coor -> 物质点坐标
% ih -> 中心点
% jh-> 邻域点
% sp-> 起始点
% ep-> 终止点
%% 作者：万冀，作者邮箱: wanji@whu.edu.cn
%% 单位：武汉大学
x = coor(:,1);
y = coor(:,2);
lineFun = @(x1,y1,x2,y2,x,y) (y-y1).*(x2-x1)-(y2-y1).*(x-x1);
xi = x(ih);
yi = y(ih);
xj = x(jh);
yj = y(jh);
fail = lineFun(sp(1),sp(2),ep(1),ep(2),xi,yi) .* ...
    lineFun(sp(1),sp(2),ep(1),ep(2),xj,yj)<=0 & ...
     lineFun(xi,yi,xj,yj,sp(1),sp(2)).*lineFun(xi,yi,xj,yj,ep(1),ep(2))<=0;
end