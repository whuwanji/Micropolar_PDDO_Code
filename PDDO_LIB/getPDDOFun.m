function gf = getPDDOFun(ndim, nSeries, coor, ih, jh, accumPos,  wfvj)
X = coor(jh,:) - coor(ih,:); % 相对位置态
q = subIndex(ndim, nSeries); % 幂指标
np = size(coor,1); % 物质点数目
b = prod(factorial(q),2); % b矩阵
A = zeros(size(q,1),size(q,1),np); %A矩阵
for i = 1:1:size(q,1)
    for j = i:1:size(q,1)
        s = wfvj;
        for k = 1:1:ndim
            s = s.*X(:,k).^(q(i,k)+q(j,k));
        end
        A(i,j,:) = accumarray2(accumPos, s);
        if(i~=j)
            A(j,i,:) = A(i,j,:);
        end
    end
end
a = zeros(size(A));
for i = 1:1:np
    a(:,:,i) = A(:,:,i)\diag(b);
end
gf = zeros(size(jh,1),size(q,1));    
unitVector = zeros(size(jh,1),size(q,1));
for j = 1:1:size(q,1)
    unitVector(:,j) = X(:,1).^(q(j,1));
    for k = 2:1:ndim
        unitVector(:,j) = unitVector(:,j).*X(:,k).^(q(j,k));
    end
end
for i = 1:1:np
  i1 = accumPos(i)+1;
  i2 = accumPos(i+1);
    for k = 1:1:size(q,1)
        gf(i1:i2,k) = unitVector(i1:i2,:)*a(:,k,i);
    end
end
end