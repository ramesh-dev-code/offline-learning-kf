%% Revised by Ramesh Perumal on Feb 17, 2016
%%
clear all;
close all;
clc;
datasets = {'R3_M1D','R3_M1U','R1_M2D','R1_M2U','R2_SD','R1_SU','R2_STRI','R1_THAL'};
t1 = 1;
t2 = 60000;
n_col = length(datasets);
x = zeros(60000,n_col);
xfn = x;
% data_loc = strcat('D:\Research\Datasets\hvs1.mat');
% tmp = importdata(data_loc);
for i = 1:n_col
    data = importdata(strcat('D:\Research\HVS\SignalProcessing\TimePointPrediction\UAAR\Code\Single_Channel_LFP_Prediction\labels\GroundTruth\val\',datasets{i},'.csv'));        
    x(:,i) = data(:,1);
end
clear tmp;

%% Preprocessing
% Parks-McClellan optimal BPF Design
fs = 1000; % Unit: Hz
fc_k = [3 5 13 15];
mags_k = [0 1 0];
devs_k = [0.001 0.001 0.001];
[order_pm, wn_pm, a_pm, win_pm] = firpmord(fc_k,mags_k,devs_k,fs);
d = firpm(order_pm, wn_pm, a_pm, win_pm);

for i = 1:n_col    
    xf = filtfilt(d,1,x(:,i)); % Filtered signal
    xfn(:,i) = MinMaxNorm(xf,0,1);     
end
n = length(xfn);
clear xf;
%% Autocorrelation function
tdmax = 1000;
temp = zeros(tdmax,n_col);
acf = temp;
pval = temp;
fzc = zeros(n_col,1);
max_in = fzc;
for i = 1:n_col
    for j = 1:tdmax
        [R,P] = corrcoef(xfn(1:n-tdmax,i),xfn(1+j:n-tdmax+j,i));
        acf(j,i) = R(1,2);
        pval(j,i) = P(1,2);
    end
    % To find First Zero crossing in autocorrelation values
    z1 = diff(sign(acf(:,i)));
    zc = find(abs(z1),2);
    fzc(i) = zc(1)+1;
    % To find the first local maximum in autocorrelation values
    [v,max_in(i)] = max(acf(100:200,i));
    max_in(i) = max_in(i) + 100;
end
opt_td = round(mean(fzc));
opt_p = round(mean(max_in));
% Large-lag standard error 
llerror = zeros(tdmax,1);
approx_var = llerror;
for j = 2:tdmax    
    approx_var(j,1) = (1/n)*(1 + 2*sum(acf(1:j-1,1).^2));
    llerror(j,1) = sqrt(approx_var(j,1));
end
    
figure('color',[1 1 1]);
plot(1:tdmax,acf,'Linewidth',2);
axis([1 tdmax -2 2])
ax = gca;
set(ax,'XTick',[1,50:50:tdmax]);
xlabel('Time lag (msec)','Fontsize',12,'Fontweight','b');
ylabel('Correlation Coefficient values','Fontsize',12,'Fontweight','b');
title('Correlogram of single-channel PD signals from different brain regions','Fontsize',14,'Fontweight','b');
line('XData', [0 tdmax], 'YData', [0 0], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','k')
t6 = text(1,-0.05,'Zero reference line');
set(t6,'Color','k','Fontsize',12,'Fontweight','b');
v1 = vline(opt_td,'b','Average First Zero Crossing');
set(v1,'Linewidth',1.5,'Linestyle','--');
line('XData', [opt_p opt_p], 'YData', [-2 2], 'LineStyle', '--', ...
    'LineWidth', 2, 'Color','b')
t6 = text(opt_p+1,-1,'Average First Local Maximum');
set(t6,'Color','b','Fontsize',12,'Fontweight','b');
% v2 = vline(opt_p,'b','Average First Local Maximum');
% set(v2,'Linewidth',1.5,'Linestyle','--');
t7 = text(20,1.7,sprintf('The optimal time delay for phase space reconstruction given by average first zero crossing is tau = %d ms', opt_td));
set(t7,'Color','k','Fontsize',14,'Fontweight','b');
t7 = text(150,1.5,sprintf('The optimal AR model order given by the average first local maximum is p = %d', opt_p));
set(t7,'Color','k','Fontsize',14,'Fontweight','b');

figure('color',[1 1 1]);
% subplot(2,2,1);
gray_color_gradients = repmat(linspace(0.9,0.4,n_col).',1,3);
axes('ColorOrder',gray_color_gradients,'NextPlot','replacechildren')
plot(1:tdmax,acf,'Linewidth',1);
line('XData', [0 tdmax], 'YData', [0 0], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','k')
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'b';
ax.XTick = [1,opt_p,500,1000];
ax.YTick = [-1,0,1];
v1 = vline(opt_p,'k','Average First Local Maximum');
set(v1,'Linewidth',1.5,'Linestyle','--');
xlabel('Time lag (ms)','Fontsize',16,'Fontweight','b');
ylabel('Autocorrelation','Fontsize',16,'Fontweight','b');
%% Partial autocorrelation function
pmax = 10; 
test = xfn;
n_test = length(test);
for i = 1:n_col
    [arcoefs,E,K] = aryule(test(:,i),pmax);
    pacf(:,i) = -K;    
end
figure('color',[1 1 1]);
p2 = subplot(2,2,1);
ax2 = gca;
plot(1:pmax,pacf,'-s');
ylim([-1,1]);
xlabel('Time lag (ms)','Fontsize',14,'Fontweight','b');
ylabel('Partial Autocorrelation','Fontsize',14,'Fontweight','b');
% xlim([0 160]);
title('Partial Autocorrelation Method','Fontsize',14,'Fontweight','b');
% 95% confidence bounds
uconf = 1.96/sqrt(n_test); 
lconf = -uconf;
axes(p2)
h2 = hline(uconf,'k');
% set(h2,'Linewidth',1.5);
h3 = hline(lconf,'k');
% set(h3,'Linewidth',1.5);
ax2.FontSize = 16;
ax2.FontWeight = 'bold';
ax2.XTick = [1:pmax];

figure('color',[1 1 1]);
% subplot(2,2,1);
gray_color_gradients = repmat(linspace(0.9,0.4,n_col).',1,3);
axes('ColorOrder',gray_color_gradients,'NextPlot','replacechildren')
p3 = plot(1:pmax,pacf,'-s');
% 95% confidence bounds
uconf = 1.96/sqrt(n_test); 
lconf = -uconf;
hold on;
h2 = hline(uconf,'k');
h3 = hline(lconf,'k');
xlim([1,pmax]);
ylim([-1,1]);
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'b';
ax.XTick = 1:pmax;
ax.YTick = [-1,-0.5,0,0.5,1];
xlabel('Lag','Fontsize',16,'Fontweight','b');
ylabel('Partial Autocorrelation','Fontsize',16,'Fontweight','b');

