function status = unitTest_factorizeRegressor

gt_r1 = [1; 3; 1];
gt_r2 = [1; 2; 1];

r = [1 5 2 6 5 1];
[r1, r2] = factorizeRegressor(r);

if norm(r1-gt_r1)+norm(r2-gt_r2)<1e-4 || norm(r2-gt_r1)+norm(r1-gt_r2)<1e-4
    status = true;
else
    status = false;
    error('factorizeRegressor is not working correctly.\n');
end

end