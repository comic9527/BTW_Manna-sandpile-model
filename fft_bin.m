function [bin_x,P_average] = fft_bin(x,P,base)

max_x = max(x);
qmin = min(x);
bin_edges = qmin*base.^(0:ceil(log(max_x/qmin)/log(base)));
segment_info = discretize(x, bin_edges);
P_average = zeros(1, length(bin_edges)-1);
for i = 1:length(bin_edges)-1
    idx = segment_info == i; 
    P_segment = P(idx); 
    P_average(i) = mean(P_segment); 
end
bin_x = bin_edges(2:end);
end
