function [data,gt,opt] = getMSRdata
    jointInd = 13;
    load ../expData/msr.mat;
    jointTraj = msr(3*jointInd-2:3*jointInd,:);
    y = jointTraj(2,:);
    y1 = y(1:16);
    [y1h,eta,v,R] = fast_incremental_hstln_mo(y1, 0.01);
    y2 = y(16:30);
    [y2h,eta,v,R] = fast_incremental_hstln_mo(y2, 0.05);
    y3 = y(30:end);
    [y3h,eta,v,R] = fast_incremental_hstln_mo(y3, 0.1);
    yn = [y2h(7:end) y3h(2:end)];
%     yn = [y1h y2h(2:end) y3h(2:end)];
    gt = [ones(1,10), 2*ones(1,22)];
    data = yn;

end