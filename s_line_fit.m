function [slope] = s_line_fit(x,y,cr)
ftx = log(x);
fty = log(y);
[xData, yData] = prepareCurveData( ftx, fty );
ft = fittype( 'poly1' );
index = floor(4/5*(length(xData)));            
[fitresult, ~] = fit( xData(1:index), yData(1:index), ft );
hold on;
scatter( xData, yData ,40,cr, 'filled' );
plot( fitresult,cr);
slope = fitresult.p1;
end