function [slope] = line_fit(x,y,cr)
ftx = log(x);
fty = log(y);
[xData, yData] = prepareCurveData( ftx, fty );
ft = fittype( 'poly1' );
[fitresult, ~] = fit( xData, yData, ft );
hold on;
scatter( xData, yData ,40,cr, 'filled' );
plot( fitresult,cr);
slope = fitresult.p1;
end