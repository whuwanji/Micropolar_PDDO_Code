function f = accumarray2(accumPos, state)
np = size(accumPos,1) - 1;
ndim = size(state, 2);
f = zeros(np, ndim);
for i = 1:1:np
    i1 = accumPos(i)+1;
    i2 = accumPos(i+1);
    f(i,:) = sum(state(i1:i2,:),1);
end
end