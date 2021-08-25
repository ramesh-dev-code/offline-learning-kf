% Programmed by Ramesh Perumal on March 27, 2019
% This script is used to predict the LFPs from the 144-sample-long LFP
% segment to increase the segment length to 512 samples. Then, 512-point
% FFT is used to compute the power spectral density and the relative HVS
% power. The LFP segment is classified as HVS if its relative HVS power 
% exceeds a predefined threshold. 
%% Online Prediction phase  
clear all;    
close all;
clc;
dataset_name = 'R1_M1U';
data = importdata(char(strcat(dataset_name,'.csv')));
% Min-Max Normalization
xn = MinMaxNorm(data(:,1),0,1);
fs = 1000; % Unit: Hz
test = xn;
N = length(test);

%% Activate this segment to load the parameters of AR Model at Interval  
p = 6; % Model Order
T = 24; % Modeling Interval
wlen = p*T; % Effective model order
slide = 24; % Unit: ms
% Trained AR_T model coefficients
Wc = [1.2834 -1.5510 1.5003 -1.1288 0.9320 -0.0405];
hvs_power = zeros(length(wlen:slide:N),1);
rel_hvs_bp = hvs_power;
%% Prediction Logic
pwin = 512; % prediction window
hamm_win = hamming(pwin);
seg = 0;
temp = zeros(pwin,1);
tic;
for ts = wlen:slide:N
    pred = temp;
    pred(1:wlen,1) = test(ts-wlen+1:ts,1);
    for k = wlen+1:1:pwin
        pred(k,1) = Wc*pred(k-T:-T:k-wlen,1);         
    end    
    % Estimation of PSDs of Predicted LFPs
    seg = seg + 1;
    [f_fft,pred_psd] = fourierpsd(hamm_win.*pred,fs);
    hvs_power(seg) = sum(pred_psd(4:7));
    % Relative HVS Power
    rel_hvs_bp(seg) = hvs_power(seg)/sum(pred_psd);
end
toc;

%% Detection Logic
% Pre-determined threshold for each LFP dataset by evaluating the detection 
% perfromance on the 1-min-long LFP dataset
tr = 0.01; 
hvs = rel_hvs_bp>tr;
t_st = wlen;
t_end = N;

% Display power features
figure('color',[1 1 1]);
subplot(2,1,1);
plot(t_st:t_end,xn(t_st:t_end),'r','Linewidth',1.5);
hold on;
plot(t_st:slide:t_end,hvs,'k','Linewidth',1.5);
xlim([1,t_end]);
ylim([-0.2,1.2]);
ylabel('Normalized Amplitude','FontSize',10,'Fontweight','b');
xlabel('Time (ms)','FontSize',10,'Fontweight','b');
title('A','FontSize',12,'Fontweight','b');

subplot(2,1,2);
plot(t_st:slide:t_end,rel_hvs_bp,'k','Linewidth',1.5);
xlim([1,t_end]);
ylabel('Relative HVS Power (%)','FontSize',10,'Fontweight','b');
xlabel('Time (ms)','FontSize',10,'Fontweight','b');
title('B','FontSize',12,'Fontweight','b');
