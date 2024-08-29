function [slope] = w_line_fit(x,y,cr)
ftx = log(x);
fty = log(y);
[xData, yData] = prepareCurveData( ftx, fty );
index = find(xData>=1/2*max(xData),1);
ft = fittype( 'poly1' );
[fitresult, ~] = fit( xData(3:index), yData(3:index), ft );
hold on;
scatter( xData, yData ,40,cr, 'filled' );
plot( fitresult,cr);
slope = fitresult.p1;
legend(['Slope: ', num2str(slope)])
end