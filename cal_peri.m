function perimeter = cal_peri(A)
boundaries = bwboundaries(A);
boundary = boundaries{1};
perimeter = 0;
for i = 1:size(boundary, 1)-1
    dx = boundary(i+1, 1) - boundary(i, 1);
    dy = boundary(i+1, 2) - boundary(i, 2);
    perimeter = perimeter + sqrt(dx^2 + dy^2);
end
dx = boundary(1, 1) - boundary(end, 1);
dy = boundary(1, 2) - boundary(end, 2);
perimeter = perimeter + sqrt(dx^2 + dy^2);
end

