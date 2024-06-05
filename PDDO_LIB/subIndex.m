function p = subIndex(M,N)
% Given an integer N, the summation of these M non-negative integers be less or equal than N.
% List all cases of p1,p2,...pM, such that sum(p1...pM) <= N.
s = (N+1)*ones(1,M);
r = (1:prod(s))';
f = cell(1,M);
[f{:}] = ind2sub(s,r);
f = [f{:}];
p = f-1;
q = sum(p,2);
[~,idx] = sort(q);
p = p(idx,:);
q = sum(p,2)<=N&sum(p,2)>=1;
p = p(q,:);
end