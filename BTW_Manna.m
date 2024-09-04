clear;clc;
N = 128;
num = 4e4;
path0 = pwd;
base = 1.1;
pro = [0 1];
for k = 1:length(pro)
    close all;
    p = pro(k); % p=0:BTW,p = 1:manna
    str = num2str(p);
    str = strrep(str,'.','p');
    path1 = fullfile(path0,['p_',str],num2str(N));
    mkdir(path1);
    sand = randi([0, 3], N);
    s = zeros(1,num);
    t = zeros(1,num);
    area = zeros(1,num);
    Rg = zeros(1,num);
    perim = zeros(1,num);
    ns = zeros(1,num);
    san3per = zeros(1,num);
    be = zeros(1,num);
    T = 0;
    clet1 = 10; small_ratio = 0.01;
    clet2 = 1;
    big_ratio = -0.35*(p-2);
    clet_matrixs1 = true; clet_matrixs2 = true;
    count =0;
    matrixs = cell(1,2);
    fg1 = zeros(1,clet1+clet2);
    jug_matrix = zeros(N);
    for n = 1:num

        si = 0;
        ti = 0;
        bei = 0;
        ar_m = zeros(N,N);
        in = unidrnd(N,1,2);
        ar_m(in(1),in(2)) = ar_m(in(1),in(2)) + 1;
        sand(in(1),in(2)) = sand(in(1),in(2)) + 1;
        new_sand = sand;
        while any(new_sand(:) >= 4)
            st = 0;
            Be(T+1) = 0;
            for i = 1:N
                for j = 1:N
                    if new_sand(i, j) >= 4
                        r = 2*randi([0,1])-1;
                        if rand <= p
                            if i-1 >= 1
                                new_sand(i-1, j) = new_sand(i-1, j) + 1 + r;
                                ar_m(i-1,j) = (1 + r)/2;
                            else
                                Be(T+1) = Be(T+1)+1+r;
                                bei = bei + 1+r;
                            end
                            if i+1 <= N
                                new_sand(i+1, j) = new_sand(i+1, j) + 1 + r;
                                ar_m(i+1,j) = (1 + r)/2;
                            else
                                Be(T+1) = Be(T+1)+1+r;
                                bei = bei + 1+r;
                            end
                            if j-1 >= 1
                                new_sand(i, j-1) = new_sand(i, j-1) + 1 - r;
                                ar_m(i,j-1) = (1 - r)/2;
                            else
                                Be(T+1) = Be(T+1)+1-r;
                                bei = bei + 1-r;
                            end
                            if j+1 <= N
                                new_sand(i, j+1) = new_sand(i, j+1) + 1 - r;
                                ar_m(i,j+1) = (1 - r)/2;
                            else
                                Be(T+1) = Be(T+1)+1-r;
                                bei = bei + 1-r;
                            end
                            new_sand(i, j) = new_sand(i, j) - 4;
                            si = si + 1;
                            st = st + 1;
                        else
                            if i-1 >= 1
                                new_sand(i-1, j) = new_sand(i-1, j) + 1;
                                ar_m(i-1,j) = 1;
                            else
                                Be(T+1) = Be(T+1)+1;
                                bei = bei+1;
                            end
                            if i+1 <= N
                                new_sand(i+1, j) = new_sand(i+1, j) + 1;
                                ar_m(i+1,j) = 1;
                            else
                                Be(T+1) = Be(T+1)+1;
                                bei = bei+1;
                            end
                            if j-1 >= 1
                                new_sand(i, j-1) = new_sand(i, j-1) + 1;
                                ar_m(i,j-1) = 1;
                            else
                                Be(T+1) = Be(T+1)+1;
                                bei = bei+1;
                            end
                            if j+1 <= N
                                new_sand(i, j+1) = new_sand(i, j+1) + 1;
                                ar_m(i,j+1) = 1;
                            else
                                Be(T+1) = Be(T+1)+1;
                                bei = bei+1;
                            end
                            new_sand(i, j) = new_sand(i, j) - 4;
                            si = si + 1;
                            st = st + 1;
                        end
                    end
                end
            end
            ti = ti + 1;
            T = T + 1; S(T) = st;
            Aves(T) = sum(new_sand(:))/(N*N);
            San3(T) = sum(new_sand(:)==3)./(N*N);
        end
        result_matrix = ar_m .* jug_matrix;
        if sum(ar_m(:))/(N*N)>=small_ratio && clet_matrixs1 && sum(result_matrix(:)) == 0
            jug_matrix = jug_matrix + ar_m;
            count = count + 1;
            fg1(count) = n;
            if count == clet1
                clet_matrixs1 = false;
            end
        end

        if sum(ar_m(:))/(N*N)>=big_ratio && clet_matrixs2
            count = count + 1;
            fg1(count) = n;
            matrixs{2} = ar_m;
            if count == clet1 + clet2
                clet_matrixs2 = false;
            end
        end

        area(n) = sum(ar_m(:));
        ns(n) = sum(new_sand(:));
        san3per(n) = sum(new_sand(:)==3)./(N*N);
        [row, col] = find(ar_m ~= 0);
        cc = bwconncomp(ar_m);
        stats = regionprops(cc,"Centroid");
        centroid =stats.Centroid;
        rc2 = ((row) - centroid(2)).^2 + ((col) - centroid(1)).^2;
        distan = sqrt(mean(rc2));
        perimeter = cal_peri(ar_m);

        figure(1),set(figure(1), 'Position', [100, 100, 800, 800]);imagesc(ar_m);
        titlestr = ['the step is ',num2str(n),newline,' R_{g} = ',num2str(distan),' loop length = ',num2str(perimeter)];
        title(titlestr,"FontSize",18);
        axis square;
        t(n) = ti;
        s(n) = si;
        Rg(n) = distan;
        perim(n) = perimeter;
        be(n) = bei;
        sand = new_sand;
        figure(2);set(figure(2), 'Position', [950, 100, 800, 800]);
        imagesc(sand);
        colorbar;
        titlestr = ['the step is ',num2str(n),newline,' t = ',num2str(ti)];
        title(titlestr,FontSize=18);
        axis square;
        pause(0.001);
    end
    matrixs{1} = jug_matrix;
    avesand = ns./(N*N);
    seg = 100;
    seg_point = num/seg;
    new_list = avesand(seg_point:seg_point:end);
    new_list = [avesand(1) new_list];
    diff_avesand = diff(new_list);
    index = find(diff_avesand < 0.01, 3, 'first');
    start_point = index(end)*seg_point;

    figure(3),set(figure(3),'position',[100,100,600,600]);
    axes('Units', 'normalized', 'Position', [.001 .001 .998 .998], 'XTickLabel', '', 'YTickLabel', '');
    imshow(cell2mat(matrixs(end)));axis square;axis off;


    figure(4),set(figure(4),'position',[100,100,600,600]);
    axes('Units', 'normalized', 'Position', [.001 .001 .998 .998], 'XTickLabel', '', 'YTickLabel', '');
    imshow(cell2mat(matrixs(1)));axis square;axis off;

    start_point2 = sum(t(1:start_point-1));
    Y1 = fft(Aves(start_point2:end));
    Q1 = abs(Y1).^2/(T-start_point2+1);
    Y2 = fft(San3(start_point2:end));
    Q2 = abs(Y2).^2/(T-start_point2+1);
    Y3 = fft(S(start_point2:end));
    Q3 = abs(Y3).^2/(T-start_point2+1);
    Y4 = fft(Be(start_point2:end));
    Q4 = abs(Y4).^2/(T-start_point2+1);
    len = length(Y1);
    a=1;
    qmin=(2.*pi)./(len*a);
    qmax=pi./a;
    q_x=linspace(qmin,qmax,len);
    [qx1,qy1] = fft_bin(q_x,Q1,base);
    [qx2,qy2] = fft_bin(q_x,Q2,base);
    [qx3,qy3] = fft_bin(q_x,Q3,base);
    [qx4,qy4] = fft_bin(q_x,Q4,base);
    figure(5),set(figure(5), 'Position', [100, 10, 2000, 800]);
    md5 = 1;
    md52 = 2;

    subplot(2,4,1);plot(1:T,Aves,'linewidth',2);xlabel('Number of time step');axis square;
    title('Average number of sand per cell');
    axchange(md5);
    subplot(2,4,5);
    plot(qx1,qy1,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md52);
    subplot(2,4,2);plot(1:T,San3,'linewidth',2);xlabel('Number of time step');axis square;
    title('Proportion of grid points with sand number 3');
    axchange(md5);
    subplot(2,4,6);plot(qx2,qy2,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md52);
    subplot(2,4,3);plot(1:T,S,'linewidth',2);xlabel('Number of time step');axis square;
    title('Number of reactive cells');
    axchange(md5);
    subplot(2,4,7);plot(qx3,qy3,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md52);
    subplot(2,4,4);plot(1:T,Be,'linewidth',2);xlabel('Number of time step');axis square;
    title('Number of sand that escapes from the boundary');
    axchange(md5);
    subplot(2,4,8);plot(qx4,qy4,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md52);

    qpsdata1 = [q_x',Q1',Q2',Q3',Q4'];
    Ta1 = array2table(qpsdata1,'VariableNames',{'qps_x','ps_Aves','ps_San3','ps_S','ps_Be'});
    T_psname1 = ('Power spectrum with the number of time step.txt');
    writetable(Ta1, fullfile(path1,T_psname1), 'Delimiter', '\t');
    qpsdata2 = [qx1',qy1',qy2',qy3',qy4'];
    Ta2 = array2table(qpsdata2,'VariableNames',{'qpx','ps_Aves_bin','ps_San3_bin','ps_S_bin','ps_Be_bin'});
    T_psname1 = ('Power spectrum with the number of time step_binning.txt');
    writetable(Ta2, fullfile(path1,T_psname1), 'Delimiter', '\t');

    flidata0 = {(1:T)',Aves',San3',S',Be'};
    max_length = max(cellfun(@numel, flidata0));
    flidata_padded0 = cellfun(@(x) [x; nan(max_length-length(x),1)], flidata0, 'UniformOutput', false);
    flidata01 = cell2mat(flidata_padded0);
    Ta3 = array2table(flidata01, 'VariableNames', {'T','Aves','San3','S','Be'});
    fliname = ('Number of time step_Steady state verification curve.txt');
    writetable(Ta3, fullfile(path1,fliname), 'Delimiter', '\t');

    Rg2 = round(Rg);
    perim2 = round(perim);
    figure(6),set(figure(6), 'Position', [100, 100, 1600, 800]);
    md6 = 2;

    ws = s(start_point:end);
    wt = t(start_point:end);
    warea = area(start_point:end);
    wRg2 = Rg2(start_point:end);
    wperim2 = perim2(start_point:end);

    ave_s = sum(ws)/(num-start_point+1);
    ave_t = sum(wt)/(num-start_point+1);
    ave_a = sum(warea)/(num-start_point+1);
    ave_Rg = sum(wRg2)/(num-start_point+1);
    ave_peri = sum(wperim2)/(num-start_point+1);

    sgtitle(['The starting point for interception is: ',num2str(start_point)]);
    subplot(2,3,1);[wsx,wsy] = whisto(ws);xlabel('s');ylabel('Proportion');axis square;axchange(md6);
    subplot(2,3,2);[wtx,wty] = whisto(wt);xlabel('t');ylabel('Proportion');axis square;axchange(md6);
    subplot(2,3,3);[wax,way] = whisto(warea);xlabel('area');ylabel('Proportion');axis square;axchange(md6);
    subplot(2,3,4);[wrx,wry] = whisto(wRg2);xlabel('R_{g}');ylabel('Proportion');axis square;axchange(md6);
    subplot(2,3,5);[wpx,wpy] = whisto(wperim2);xlabel('peri');ylabel('Proportion');axis square;axchange(md6);


    wsy_acc = cumsum(wsy(end:-1:1));
    wsy_acc = wsy_acc(end:-1:1);
    wty_acc = cumsum(wty(end:-1:1));
    wty_acc = wty_acc(end:-1:1);
    way_acc = cumsum(way(end:-1:1));
    way_acc = way_acc(end:-1:1);
    wry_acc = cumsum(wry(end:-1:1));
    wry_acc = wry_acc(end:-1:1);
    wpy_acc = cumsum(wpy(end:-1:1));
    wpy_acc = wpy_acc(end:-1:1);
    figure(7);set(figure(7), 'Position', [100, 100, 1600, 800]);
    md7 = 2;
    subplot(2,3,1);plot(wsx,wsy_acc,'linewidth',2);xlabel('s');ylabel('Cumulative frequency');axis square;axchange(md7);
    subplot(2,3,2);plot(wtx,wty_acc,'linewidth',2);xlabel('t');ylabel('Cumulative frequency');axis square;axchange(md7);
    subplot(2,3,3);plot(wax,way_acc,'linewidth',2);xlabel('area');ylabel('Cumulative frequency');axis square;axchange(md7);
    subplot(2,3,4);plot(wrx,wry_acc,'linewidth',2);xlabel('R_{g}');ylabel('Cumulative frequency');axis square;axchange(md7);
    subplot(2,3,5);plot(wpx,wpy_acc,'linewidth',2);xlabel('peri');ylabel('Cumulative frequency');axis square;axchange(md7);


    wdata = {wsx',wsy',wtx',wty',wax',way',wrx',wry',wpx',wpy'};
    adata = {wsx',wsy_acc',wtx',wty_acc',wax',way_acc',wrx',wry_acc',wpx',wpy_acc'};
    max_length = max(cellfun(@numel, wdata));
    max_length2 = max(cellfun(@numel, adata));
    wdata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata, 'UniformOutput', false);
    wdata1 = cell2mat(wdata_padded);
    adata_padded = cellfun(@(x) [x; nan(max_length2-length(x),1)], adata, 'UniformOutput', false);
    adata1 = cell2mat(adata_padded);

    Ta4 = array2table(wdata1, 'VariableNames', {'wsx','wsy','wtx','wty','wax','way','wrx','wry','wpx','wpy'});
    wname = ('Frequency curve.txt');
    writetable(Ta4, fullfile(path1,wname), 'Delimiter', '\t');
    Ta4_a = array2table(adata1, 'VariableNames', {'wsx','wsy_acc','wtx','wty_acc','wax','way_acc','wrx','wry_acc','wpx','wpy_acc'});
    aname = [num2str(N), 'Critical height accumulation.txt'];
    writetable(Ta4_a, fullfile(path1,aname), 'Delimiter', '\t');



    figure(8),set(figure(8), 'Position', [100, 100, 1600, 800]);
    md8 = 1;
    subplot(2,3,1);wssl= w_line_fit(wsx,wsy,'b');xlabel('log(s)');ylabel('log(Frequency)');axis square;hold off;axchange(md8);
    subplot(2,3,2);wtsl = w_line_fit(wtx,wty,'b');xlabel('log(t)');ylabel('log(Frequency)');axis square;hold off;axchange(md8);
    subplot(2,3,3);wasl = w_line_fit(wax,way,'b');xlabel('log(area)');ylabel('log(Frequency)');axis square;hold off;axchange(md8);
    subplot(2,3,4);wrsl = w_line_fit(wrx,wry,'b');xlabel('log R_{g}', 'Interpreter', 'tex');ylabel('log(Frequency)');axis square;hold off;axchange(md8);
    subplot(2,3,5);wpsl = w_line_fit(wpx,wpy,'b');xlabel('log(peri)');ylabel('log(Frequency)');axis square;hold off;axchange(md8);

    figure(9),set(figure(9), 'Position', [100, 10, 2000, 800]);
    md10 = 1;
    md92 = 2;
    X1 = fft(avesand(start_point:end));
    P1 = abs(X1).^2/(num-start_point+1);
    X2 = fft(san3per(start_point:end));
    P2 = abs(X2).^2/(num-start_point+1);
    X3 = fft(s(start_point:end));
    P3 = abs(X3).^2/(num-start_point+1);
    X4 = fft(be(start_point:end));
    P4 = abs(X4).^2/(num-start_point+1);
    len1 = length(X1);
    a1=1;
    qmin1=(2.*pi)./(len1*a1);
    qmax1=pi./a1;
    p_x=linspace(qmin1,qmax1,len1);
    [px1,py1] = fft_bin(p_x,P1,base);
    [px2,py2] = fft_bin(p_x,P2,base);
    [px3,py3] = fft_bin(p_x,P3,base);
    [px4,py4] = fft_bin(p_x,P4,base);

    subplot(2,4,1);plot(1:num,avesand,'linewidth',2);xlabel('Number of input sand');axis square;
    title('Average number of sand per cell');
    axchange(md10);
    subplot(2,4,5);plot(px1,py1,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md92);
    subplot(2,4,2);
    plot(1:num,san3per,'linewidth',2);xlabel('Number of input sand');axis square;
    title('Proportion of grid points with sand number 3');
    axchange(md10);
    subplot(2,4,6);plot(px2,py2,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md92);
    subplot(2,4,3);plot(1:num,s,'linewidth',2);xlabel('Number of input sand');axis square;
    title('Number of reactive cells');
    axchange(md10);
    subplot(2,4,7);plot(px3,py3,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md92);
    subplot(2,4,4);plot(1:num,be,'linewidth',2);xlabel('Number of input sand');axis square;
    title('Number of sand that escapes from the boundary');
    axchange(md10);
    subplot(2,4,8);
    plot(px4,py4,'r','linewidth',2);hold off;axis square;
    title('Power spectrum');xlabel('Frequency');ylabel('Intensity');
    axchange(md92);

    ppsdata1 = [p_x',P1',P2',P3',P4'];
    Ta5 = array2table(ppsdata1,'VariableNames',{'pps_x','ps_avesand','ps_san3per','ps_s','ps_be'});
    N_psname1 = ('Power spectrum with the addition of sand particles.txt');
    writetable(Ta5, fullfile(path1,N_psname1), 'Delimiter', '\t');
    ppsdata2 = [px1',py1',py2',py3',py4'];
    Ta6 = array2table(ppsdata2,'VariableNames',{'ppx','ps_avesand_bin','ps_san3per_bin','ps_s_bin','ps_be_bin'});
    T_psname1 = ('Power spectrum with the addition of sand particles_binning.txt');
    writetable(Ta6, fullfile(path1,T_psname1), 'Delimiter', '\t');

    flidata = {(1:num)',avesand',san3per',s',be'};
    max_length = max(cellfun(@numel, flidata));
    flidata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], flidata, 'UniformOutput', false);
    flidata1 = cell2mat(flidata_padded);
    Ta7 = array2table(flidata1, 'VariableNames', {'x','avesand','san3per','s','be'});
    fliname = ('Number of sand particles_Steady state verification curve.txt');
    writetable(Ta7, fullfile(path1,fliname), 'Delimiter', '\t');

    wRg = Rg(start_point:end);
    wperim = perim(start_point:end);
    [Rg9, idx] = sort(wRg);
    s9 = ws(idx);
    t9 = wt(idx);
    area9 = warea(idx);
    perim9 = wperim(idx);

    md10 = 2;
    [unique_R, ~, idx_unique] = unique(Rg9);
    Rg10 = unique_R;
    s10 = accumarray(idx_unique, s9, [], @mean);
    t10 = accumarray(idx_unique, t9, [], @mean);
    area10 = accumarray(idx_unique, area9, [], @mean);
    perim10 = accumarray(idx_unique, perim9, [], @mean);
    figure(10);sgtitle('Scaling relations');
    set(figure(10), 'Position', [100, 100, 1200, 600]);
    [~,s_average]= binning(Rg10,s10,base);
    [~,t_average]=binning(Rg10,t10,base);
    [~,area_average]=binning(Rg10,area10,base);
    [bins,perim_average]=binning(Rg10,perim10,base);

    subplot(1,2,1);
    plot(bins,s_average,'-o',bins,t_average,'-o',bins,area_average,'-o',bins,perim_average,'-o','LineWidth',2);
    legend('s','t','area','peri','Location', 'NorthWest');xlabel('R_{g}');axis square;
    axchange(md10);

    subplot(1,2,2);
    sl1 = s_line_fit(bins,s_average,'c');
    sl2 = s_line_fit(bins,t_average,'b');
    sl3 = s_line_fit(bins,area_average,'r');
    sl4 = s_line_fit(bins,perim_average,'black');
    hold off;
    legend( {'average s', ['Fitting of average s: ','slope = ',num2str(sl1)], 'average t', ['Fitting of average t: ','slope = ',num2str(sl2)],'average area', ['Fitting of average area: ','slope = ',num2str(sl3)],'average perim', ['Fitting of average perim: ','slope = ',num2str(sl4)]} ...
        , 'Location', 'NorthWest', 'Interpreter', 'none' );
    xlabel('log R_{g}', 'Interpreter', 'tex');
    ylabel( 'log(y)', 'Interpreter', 'none' )
    axis square;

    ssandata = {Rg10',s10,t10,area10,perim10,log(bins'),log(s_average'),log(t_average'),log(area_average'),log(perim_average')};
    max_length = max(cellfun(@numel, ssandata));
    ssandata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], ssandata, 'UniformOutput', false);
    ssandata1 = cell2mat(ssandata_padded);
    Ta8 = array2table(ssandata1, 'VariableNames', {'Rgo','so','to','areao','perimo','Rg','s','t','area','perim'});
    ssanname = ('Scaling relations.txt');
    writetable(Ta8, fullfile(path1,ssanname), 'Delimiter', '\t');

    table_single_vb = table({'ave_s';'ave_t';'ave_a';'ave_Rg';'ave_peri';'start_point_N'; 'start_point2_T';'T'; 'base'}, [ave_s;ave_t;ave_a;ave_Rg;ave_peri;start_point;start_point2;T;base], 'VariableNames', {'Label', 'Value'});
    writetable(table_single_vb, fullfile(path1,'table_single_vb.txt'), 'Delimiter', '\t');
    ordata1 = [(1:num)',s',t',area',Rg',perim',ns',san3per',be'];
    ordata2 = [(1:T)',Aves',San3',S',Be'];
    Ta9 = array2table(ordata1, 'VariableNames', {'num','s','t','area','Rg','perim','ns','san3per','be'});
    Ta10 = array2table(ordata2, 'VariableNames', {'T','Aves','San3','S','Be'});
    or1name = ('origindata_N.txt');
    or2name = ('origindata_T.txt');
    writetable(Ta9, fullfile(path1,or1name), 'Delimiter', '\t');
    writetable(Ta10, fullfile(path1,or2name), 'Delimiter', '\t');

    sldata1 = [wssl,wtsl,wasl,wrsl,wpsl,sl1,sl2,sl3,sl4];
    Ta11 = array2table(sldata1,'VariableNames',{'wssl','wtsl','wasl','wrsl','wpsl','s_dispersion slope','t_dispersion slope','area_dispersion slope','perim_dispersion slope'});
    slname = ('Frequency slope and dispersion slope.txt');
    writetable(Ta11, fullfile(path1,slname), 'Delimiter', '\t');

    allFigures = get(0, 'Children');
    numFiguresToSave = length(allFigures);
    for i = 1:numFiguresToSave
        figHandle = allFigures(i);
        figPosition = get(figHandle, 'Position');
        num1 = 4;
        num2 = num1*figPosition(3)/figPosition(4);
        filename = sprintf('figure%d.png',numFiguresToSave-i+1);
        saveas(figHandle, fullfile(path1,filename), 'png');
    end
    clear ('Aves','San3','Be','S');
end
