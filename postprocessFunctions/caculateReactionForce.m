function [u1, u2, fx, fy] = caculateReactionForce(folderName,startStep, step, endStep)
a=load([folderName,'\pd_model.mat']);
coor = a.pd_model.coor;
stepss = startStep:step:endStep;
fx = zeros(size(stepss));
fy = zeros(size(fx));
u1 = zeros(size(fx));
u2 = zeros(size(fx));
downNode = find(coor(:,2)<=-0.5+eps);
ct=1;
for i = startStep:step:endStep
    ct = ct+1;
    resFile = [folderName,'\result_step_',num2str(i),'.mat'];
    if(exist(resFile,'file'))
        b = load(resFile);
    else
        break;
    end
    s21 = full(b.result.s21);
    s22 = full(b.result.s22);
    dx = 1/101;
    fx(ct) = sum(-s21(downNode)*dx*dx);
    fy(ct) = sum(-s22(downNode)*dx*dx);
    u1(ct) = b.result.apply_disp_u1;
    u2(ct) = b.result.apply_disp_u2;
end

