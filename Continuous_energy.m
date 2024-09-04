close all;clear;clc;
tic;
N = 32;
num =8000;
path0=pwd;
path1 = fullfile(path0,num2str(N));
mkdir(path1);

energy = zeros(N);
E_matrix=zeros(num,N,N);

rs=round(rand(1,num) * 100) / 100;
s = zeros(1,num);
t = zeros(1,num);
area = zeros(1,num);
Rg = zeros(1,num);
perim = zeros(1,num);
ns = zeros(1,num);
be = zeros(1,num);
T = 0;
base=1.1;

clet1 = 10;
small_ratio = 0.01;
clet2 = 1;
big_ratio = 0.7;
clet_matrixs1 = true; clet_matrixs2 = true;
count =0;
matrixs = cell(1,2);
fg1 = zeros(1,clet1+clet2);
jug_matrix=zeros(N);


for n = 1:num
    n

    si = 0;
    ti = 0;
    bei = 0;

    ar_m = zeros(N,N);
    in = unidrnd(N,1,2);
    ar_m(in(1),in(2)) = ar_m(in(1),in(2)) + 1;
    energy(in(1),in(2)) = energy(in(1),in(2)) + rs(n);
    new_energy = energy;
    while any(new_energy(:) >= 1)
        st = 0;
        Be(T+1) = 0;

        for i = 1:N
            for j = 1:N

                if new_energy(i, j) >= 1

                    if i-1 >= 1
                        new_energy(i-1, j) = new_energy(i-1, j) + new_energy(i, j)/4;
                        ar_m(i-1,j) = 1;
                    else
                        Be(T+1) = Be(T+1)+new_energy(i, j)/4;
                        bei = bei + new_energy(i, j)/4;
                    end
                    if i+1 <= N
                        new_energy(i+1, j) = new_energy(i+1, j) + new_energy(i, j)/4;
                        ar_m(i+1,j) = 1;
                    else
                        Be(T+1) = Be(T+1)+new_energy(i, j)/4;
                        bei = bei + new_energy(i, j)/4;
                    end
                    if j-1 >= 1
                        new_energy(i, j-1) = new_energy(i, j-1) + new_energy(i, j)/4;
                        ar_m(i,j-1) = 1;
                    else
                        Be(T+1) = Be(T+1)+new_energy(i, j)/4;
                        bei = bei + new_energy(i, j)/4;
                    end
                    if j+1 <= N
                        new_energy(i, j+1) = new_energy(i, j+1) + new_energy(i, j)/4;
                        ar_m(i,j+1) = 1;
                    else
                        Be(T+1) = Be(T+1)+new_energy(i, j)/4;
                        bei = bei + new_energy(i, j)/4;
                    end


                    new_energy(i, j) = 0;
                    si = si + 1;
                    st = st + 1;
                end
            end
        end
        ti = ti + 1;
        T = T + 1; S(T) = st;
        avesand(T) = sum(new_energy(:))./(N*N);
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

    E_matrix(n,:,:)=new_energy;
    area(n) = sum(ar_m(:));
    ns(n) = sum(new_energy(:));

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

    energy = new_energy;
    figure(2);set(figure(2), 'Position', [950, 100, 800, 800]);
    imagesc(energy);
    colorbar;
    titlestr = ['the step is ',num2str(n),newline,' t = ',num2str(ti)];
    title(titlestr,FontSize=18);
    axis square;
    pause(0.001);
end
toc;
avesand1 = ns./(N*N);
seg = 100;
seg_point = num/seg;
new_list = avesand1(seg_point:seg_point:end);
new_list = [avesand1(1) new_list];
diff_avesand = diff(new_list);
index = find(diff_avesand < 0.003, 3, 'first');
energy_step = index(end)*seg_point;
T_start_step= sum(t(1:energy_step-1));

matrixs{1} = jug_matrix;
figure(3),set(figure(3),'position',[100,100,600,600]);
ha1 = tight_subplot(1,1,[.001 .001],[.001 .001],[.001 .001]);
imshow(cell2mat(matrixs(end)));
axis square;axis off;

figure(4),set(figure(4),'position',[100,100,600,600]);
ha2 = tight_subplot(1,1,[.001 .001],[.001 .001],[.001 .001]);
imshow(cell2mat(matrixs(1)));axis square;axis off;


[Rg_1, idx] = sort(Rg(energy_step:end));
t0=t(energy_step:end);
t1 = t0(idx);
area0=area(energy_step:end);
area1 = area0(idx);
s0=s(energy_step:end);
s1=s0(idx);
perim0=perim(energy_step:end);
perim1=perim0(idx);

[unique_Rg, ~, idx_unique] = unique(Rg_1);
Rg_2 = unique_Rg;
t2 = accumarray(idx_unique, t1, [], @mean);
area2 = accumarray(idx_unique, area1, [], @mean);
s2 = accumarray(idx_unique, s1, [], @mean);
perim2 = accumarray(idx_unique, perim1, [], @mean);
figure(5);sgtitle('Scaling relations');
set(figure(5), 'Position', [100, 100, 1200, 600]);

[~,s_average]= binning(Rg_2,s2,base);
[~,t_average]=binning(Rg_2,t2,base);
[~,area_average]=binning(Rg_2,area2,base);
[bins,perim_average]=binning(Rg_2,perim2,base);

subplot(1,2,1);
plot(bins,s_average,'-o',bins,t_average,'-o',bins,area_average,'-o',bins,perim_average,'-o','LineWidth',2);
legend('s','t','area','peri','Location', 'NorthWest');xlabel('R_{g}');axis square;
axchange(2);

subplot(1,2,2);
sl1 = line_fit(bins,s_average,'c');
sl2 = line_fit(bins,t_average,'b');
sl3 = line_fit(bins,area_average,'r');
sl4 = line_fit(bins,perim_average,'black');
hold off;

legend( {'average s', ['Fitting of average s: ','slope = ',num2str(sl1)], 'average t', ['Fitting of average t: ','slope = ',num2str(sl2)],'average area', ['Fitting of average area: ','slope = ',num2str(sl3)],'average perim', ['Fitting of average perim: ','slope = ',num2str(sl4)]} ...
    , 'Location', 'NorthWest', 'Interpreter', 'none' );
xlabel( 'log R_{g}', 'Interpreter', 'tex' );
ylabel( 'log y', 'Interpreter', 'none' );
axis square;

ssandata = {Rg_2',s2,t2,area2,perim2,log(bins'),log(s_average'),log(t_average'),log(area_average'),log(perim_average')};
max_length = max(cellfun(@numel, ssandata));
ssandata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], ssandata, 'UniformOutput', false);
ssandata1 = cell2mat(ssandata_padded);
Ta1 = array2table(ssandata1, 'VariableNames', {'Rgo','so','to','areao','perimo','Rg','s','t','area','perim'});
ssanname = ('Scaling relations.txt');
writetable(Ta1, fullfile(path1,ssanname), 'Delimiter', '\t');

Rg2 = round(Rg);
perim2=round(perim);
md5 = 1;
figure(6),set(figure(6), 'Position', [100, 100, 1600, 800]);
he_ws = sum(s(energy_step:end));
he_wt = sum(t(energy_step:end));
he_warea = sum(area(energy_step:end));
he_wRg2 = sum(Rg2(energy_step:end));
he_wperim2 = sum(perim2(energy_step:end));

ave_s=he_ws/(num-energy_step+1);
ave_T = he_wt/(num-energy_step+1);
ave_a=he_warea/(num-energy_step+1);
ave_r=he_wRg2/(num-energy_step+1);
ave_p=he_wperim2/(num-energy_step+1);

ws = s(energy_step:end);
wt = t(energy_step:end);
warea = area(energy_step:end);
wRg2 = Rg2(energy_step:end);
wperim2 = perim2(energy_step:end);
sgtitle(['The starting point for interception is ',num2str(energy_step)]);

md4 = 2;
subplot(2,3,1);[wsx,wsy]=whisto(ws);xlabel('s');ylabel('Frequency');axis square;axchange(md4);
subplot(2,3,2);[wtx,wty]=whisto(wt);xlabel('t');ylabel('Frequency');axis square;axchange(md4);
subplot(2,3,3);[wax,way]=whisto(warea);xlabel('area');ylabel('Frequency');axis square;axchange(md4);
subplot(2,3,4);[wrx,wry]=whisto(wRg2);xlabel('R_{g}');ylabel('Frequency');axis square;axchange(md4);
subplot(2,3,5);[wpx,wpy]=whisto(wperim2);xlabel('peri');ylabel('Frequency');axis square;axchange(md4);

wsx1 = wsx(~isnan(wsx));
wsy1 = wsy(~isnan(wsy));
wtx1 = wtx(~isnan(wtx));
wty1 = wty(~isnan(wty));
wax1 = wax(~isnan(wax));
way1 = way(~isnan(way));
wrx1 = wrx(~isnan(wrx));
wry1 = wry(~isnan(wry));
wpx1 = wpx(~isnan(wpx));
wpy1 = wpy(~isnan(wpy));

wsy_acc = cumsum(wsy1(end:-1:1));
wsy_acc = wsy_acc(end:-1:1);
wty_acc = cumsum(wty1(end:-1:1));
wty_acc = wty_acc(end:-1:1);
way_acc = cumsum(way1(end:-1:1));
way_acc = way_acc(end:-1:1);
wry_acc = cumsum(wry1(end:-1:1));
wry_acc = wry_acc(end:-1:1);
wpy_acc = cumsum(wpy1(end:-1:1));
wpy_acc = wpy_acc(end:-1:1);

figure(7);set(figure(7), 'Position', [100, 100, 1600, 800]);

subplot(2,3,1);plot(wsx1,wsy_acc,'linewidth',2);xlabel('s');ylabel('Cumulative frequency');axis square;axchange(2);
subplot(2,3,2);plot(wtx1,wty_acc,'linewidth',2);xlabel('t');ylabel('Cumulative frequency');axis square;axchange(2);
subplot(2,3,3);plot(wax1,way_acc,'linewidth',2);xlabel('area');ylabel('Cumulative frequency');axis square;axchange(2);
subplot(2,3,4);plot(wrx1,wry_acc,'linewidth',2);xlabel('R_{g}');ylabel('Cumulative frequency');axis square;axchange(2);
subplot(2,3,5);plot(wpx1,wpy_acc,'linewidth',2);xlabel('peri');ylabel('Cumulative frequency');axis square;axchange(2);

wdata = {wsx1',wsy_acc',wtx1',wty_acc',wax1',way_acc',wrx1',wry_acc',wpx1',wpy_acc'};
max_length = max(cellfun(@numel, wdata));
wdata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata, 'UniformOutput', false);
wdata1 = cell2mat(wdata_padded);
Ta2 = array2table(wdata1, 'VariableNames', {'wsx1','wsy_acc','wtx1','wty_acc','wax1','way_acc','wrx1','wry_acc','wpx1','wpy_acc'});
wname = 'Critical height accumulation.txt';
writetable(Ta2, fullfile(path1,wname), 'Delimiter', '\t');
wdata2 = {wsx',wsy',wtx',wty',wax',way',wrx',wry',wpx',wpy'};
max_length = max(cellfun(@numel, wdata2));
wdata_padded1 = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata2, 'UniformOutput', false);
wdata3 = cell2mat(wdata_padded1);
Ta3 = array2table(wdata1, 'VariableNames', {'wsx','wsy','wtx','wty','wax','way','wrx','wry','wpx','wpy'});
wname1 = ('Steady state verification.txt');
writetable(Ta3, fullfile(path1,wname1), 'Delimiter', '\t');


wdata4 ={ave_s',ave_T',ave_a',ave_r'};
max_length = max(cellfun(@numel, wdata4));
wdata_padded2 = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata4, 'UniformOutput', false);
wdata5 = cell2mat(wdata_padded2);
Ta4 = array2table(wdata5, 'VariableNames', {'ave_s','ave_T','ave_a','ave_r'});
wname2 =  'Average physical quantity.txt';
writetable(Ta4, fullfile(path1, wname2), 'Delimiter', '\t');

figure(8),set(figure(8), 'Position', [100, 100, 1600, 800]);
md5 = 1;
subplot(2,3,1);wssl= w_line_fit(wsx,wsy,'b');xlabel('log(s)');ylabel('log(Frequency)');axis square;hold off;axchange(md5);
subplot(2,3,2);wtsl = w_line_fit(wtx,wty,'b');xlabel('log(t)');ylabel('log(Frequency)');axis square;hold off;axchange(md5);
subplot(2,3,3);wasl = w_line_fit(wax,way,'b');xlabel('log(area)');ylabel('log(Frequency)');axis square;hold off;axchange(md5);
subplot(2,3,4);wrsl = w_line_fit(wrx,wry,'b');xlabel('log R_{g}', 'Interpreter', 'tex');ylabel('log(Frequency)');axis square;hold off;axchange(md5);
subplot(2,3,5);wpsl = w_line_fit(wpx,wpy,'b');xlabel('log(peri)');ylabel('log(Frequency)');axis square;hold off;axchange(md5);

sldata1 = [wssl,wtsl,wasl,wrsl,wpsl,sl1,sl2,sl3,sl4];
Ta5 = array2table(sldata1,'VariableNames',{'wssl','wtsl','wasl','wrsl','wpsl','sl1','sl2','sl3','sl4'});
slname = ('Frequency slope and dispersion slope.txt');
writetable(Ta5, fullfile(path1,slname), 'Delimiter', '\t');

figure(9),set(figure(9), 'Position', [100, 10, 2000, 800]);
md6 = 1;
X1 = fft(avesand1(energy_step:end));
P1 = abs(X1).^2/(num-energy_step+1);
X2 = fft(s(energy_step:end));
P2 = abs(X2).^2/(num-energy_step+1);
X3 = fft(be(energy_step:end));
P3 = abs(X3).^2/(num-energy_step+1);
p_x = 1:length(X1);

subplot(2,3,1);plot(1:num,avesand1,'linewidth',2);xlabel('Number of times input energy');axis square;
title('Average energy per site');
axchange(md6);

subplot(2,3,2);
plot(1:num,s,'linewidth',2);xlabel('Number of times input energy');axis square;
title('Active site number');
axchange(md6);

subplot(2,3,3);plot(1:num,be,'linewidth',2);xlabel('Number of times input energy');axis square;
title('Dissipated energy');
axchange(md6);

len = length(p_x);
a=1;
qmin=(2.*pi)./(len*a);
qmax=pi./a;
q=linspace(qmin,qmax,len);
x=q(1:end-1);
ux=x';
[~,pps_ave] = fft_bin(ux,P1,base,qmin);
[~,pps_s]=fft_bin(ux,P2,base,qmin);
[binn,pps_be]=fft_bin(ux,P3,base,qmin);

subplot(2,3,4);plot(binn,pps_ave,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axchange(2);
subplot(2,3,5);plot(binn,pps_s,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axchange(2);
subplot(2,3,6);plot(binn,pps_be,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axchange(2);

wdata8 = {binn',pps_ave',pps_s',pps_be'};
max_length = max(cellfun(@numel, wdata8));
wdata_padded7 = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata8, 'UniformOutput', false);
wdata9= cell2mat(wdata_padded7);
Ta16 = array2table(wdata9, 'VariableNames', {'binn','pps_ave','pps_s','pps_be'});
wname7 =  'Power spectrum with the number of energy releases_binning.txt';
writetable(Ta16, fullfile(path1, wname7), 'Delimiter', '\t');

ppsdata1 = [p_x',P1',P2',P3'];
Ta6 = array2table(ppsdata1,'VariableNames',{'pps_x','ps_avesand1','ps_s','ps_be'});
N_psname1 = ('Power spectrum with the number of energy releases.txt');
writetable(Ta6, fullfile(path1,N_psname1), 'Delimiter', '\t');

Y1 = fft(avesand(T_start_step:end));
Q1 = abs(Y1).^2/(T-T_start_step+1);
Y2 = fft(S(T_start_step:end));
Q2 = abs(Y2).^2/(T-T_start_step+1);
Y3 = fft(Be(T_start_step:end));
Q3 = abs(Y3).^2/(T-T_start_step+1);
q_x = 1:length(Y1);

figure(10),set(figure(10), 'Position', [100, 10, 2000, 800]);
md3 = 1;
subplot(2,3,1);plot(1:T,avesand,'linewidth',2);xlabel('Number of time step');
title('Average energy per site');axis square;
axchange(md3);

subplot(2,3,2);plot(1:T,S,'linewidth',2);xlabel('Number of time step');axis square;
title('Active site number');
axchange(md3);

subplot(2,3,3);plot(1:T,Be,'linewidth',2);xlabel('Number of time step');axis square;
title('Dissipated energy');
axchange(md3);

len = length(q_x);
a=1;
qmin=(2.*pi)./(len*a);
qmax=pi./a;
q=linspace(qmin,qmax,len);
x=q(1:end-1);
ux=x';
[~,ps_ave] = fft_bin(ux,Q1,base,qmin);
[~,ps_s]=fft_bin(ux,Q2,base,qmin);
[bins,ps_be]=fft_bin(ux,Q3,base,qmin);

subplot(2,3,4);plot(bins,ps_ave,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axis square;axchange(2);
subplot(2,3,5);plot(bins,ps_s,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axis square;axchange(2);
subplot(2,3,6);plot(bins,ps_be,'linewidth',2);title('Power spectrum');xlabel('Frequency');ylabel('Intensity');axis square;axchange(2);

wdata6 = {bins',ps_ave',ps_s',ps_be'};
max_length = max(cellfun(@numel, wdata6));
wdata_padded3 = cellfun(@(x) [x; nan(max_length-length(x),1)], wdata6, 'UniformOutput', false);
wdata7 = cell2mat(wdata_padded3);
Ta8 = array2table(wdata7, 'VariableNames', {'bins','ps_ave','ps_s','ps_be'});
wname3 = 'Power spectrum with the number of time step_binning.txt';
writetable(Ta8, fullfile(path1, wname3), 'Delimiter', '\t');

qpsdata1 = [q_x',Q1',Q2',Q3'];
Ta9 = array2table(qpsdata1,'VariableNames',{'qps_x','ps_Aves','ps_S','ps_Be'});
T_psname3 = ('Power spectrum with the number of time step.txt');
writetable(Ta9, fullfile(path1,T_psname3), 'Delimiter', '\t');

flidata0 = {(1:T)',avesand',S',Be'};
max_length = max(cellfun(@numel, flidata0));
flidata_padded0 = cellfun(@(x) [x; nan(max_length-length(x),1)], flidata0, 'UniformOutput', false);
flidata01 = cell2mat(flidata_padded0);
Ta11 = array2table(flidata01, 'VariableNames', {'T','Aves','S','Be'});
fliname = ('Number of time step_Steady state verification curve.txt');
writetable(Ta11, fullfile(path1,fliname), 'Delimiter', '\t');

figure(11);
CA=E_matrix(energy_step:energy_step+4000,:,:);
AA=CA(:);
his=histogram(AA, 'FaceColor', 'b', 'Normalization', 'probability');
xlabel('E');
ylabel('P(E)');
title('The energy distribution of a single node P(E) in critical states');


figure(12);
values = his.Values;
edges = his.BinEdges;
centers = (edges(1:end-1) + edges(2:end)) / 2;

plot(centers, values, 'LineWidth', 2);
xlabel('E');
ylabel('P(E)');
title('The energy distribution of a single node P(E) in critical states');

flidata = {(1:num)',avesand1',s',be'};
max_length = max(cellfun(@numel, flidata));
flidata_padded = cellfun(@(x) [x; nan(max_length-length(x),1)], flidata, 'UniformOutput', false);
flidata1 = cell2mat(flidata_padded);
Ta13 = array2table(flidata1, 'VariableNames', {'x','avesand1','s','be'});
fliname = ('Number of times input energy_Steady state verification curve.txt');
writetable(Ta13, fullfile(path1,fliname), 'Delimiter', '\t');

flidata3 = {base',T',T_start_step',energy_step',fg1',centers',values'};
max_length = max(cellfun(@numel, flidata3));
flidata_padde = cellfun(@(x) [x; nan(max_length-length(x),1)], flidata3, 'UniformOutput', false);
flidata4 = cell2mat(flidata_padde);
Ta14 = array2table(flidata4, 'VariableNames', {'base','T','T_start_step','energy_step','fg1','centers','values'});
fliname2 = ('origindata1.txt');
writetable(Ta14, fullfile(path1,fliname2), 'Delimiter', '\t');

ordata1 = [s',t',area',Rg',perim',ns',be'];
ordata2 = [S',Be'];
Ta15 = array2table(ordata1, 'VariableNames', {'s','t','area','Rg','perim','ns','be'});
Ta16 = array2table(ordata2, 'VariableNames', {'S','Be'});
or1name = ('origindata2.txt');
or2name = ('origindata3.txt');
writetable(Ta15, fullfile(path1,or1name), 'Delimiter', '\t');
writetable(Ta16, fullfile(path1,or2name), 'Delimiter', '\t');