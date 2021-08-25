%% Programmed by Ramesh Perumal
% Learning phase of the Offline-Learning Method
clear all;
close all;
clc;
% Input Raw LFP
data = importdata('R1_M1D.csv');
x = data(:,1);
N = length(x);

% Preprocessing
% Band-pass Filtering
% Parks-McClellan optimal BPF Design
fs = 1000; % Unit: Hz
fc_k = [3 5 13 15];
mags_k = [0 1 0];
devs_k = [0.001 0.001 0.001];
[order_pm, wn_pm, a_pm, win_pm] = firpmord(fc_k,mags_k,devs_k,fs);
d = firpm(order_pm, wn_pm, a_pm, win_pm);
xf = filtfilt(d,1,x); % Filtered signal
% Min-Max Normalization
xfn = MinMaxNorm(xf,0,1);

%% Initialization
p = 6; % AR model order 
T = 24; % Sampling time-lag 
wlen = p*T; % Window Length
temp1 = zeros(p,1);
X = temp1; % p x 1 AR model input vector
W = temp1; % p x 1 AR model coefficient vector 
% p x N AR model coefficient matrix to store the W values during learning
wt = zeros(p,N);    
I = eye(p); % p x p identity matrix
S = I; % State estimate error covariance matrix 
e = zeros(N,1); % Measurement residual
x_est = e; % Predicted LFPs 

% Estimation of Measurement noise variance from 1-s-long non-spindling LFP
R = var(diff(xfn(56001:57000)));
% Process noise variance optimized after several trials
cov_Q = 0.1*R*I; 
%% Kalman Filter Algorithm
% Training Phase to determine the converged W
tic;
for n = wlen+1:N
    % Prediction step of Kalman filter 
    X(:,1) = xfn(n-T:-T:n-wlen,1); % Input signal vector
    x_est(n) = W'*X;  
    S = S + cov_Q; 
    % Correction step of Kalman filter
    SX = S*X;
    temp2 = X'*SX + R;
    K = SX./temp2; % Kalman gain 
    e(n) = xfn(n) - x_est(n);     
    W = W + K*e(n); % Aposteriori estimate of W
    wt(:,n)= W; 
    S = S-(K*SX'); % % Aposteriori estimate of S
end
toc;
Wc = W; % Converged W after 60 s
%%
% Fourier spectrum of raw LFPs
[f_fft,x_psd] = fourierpsd(x,fs);
% Express the PSD in dB
x_psd_db =10*log10(x_psd);

% Frequency respose of Wc
fs_2 = 42;
b = 1;
a = [1,-1*Wc'];
N2 = 256;
[h,f] = freqz(b,a,N2,fs_2);
h_mag = abs(h).^2;
h_mag_db = 10*log10(h_mag);

% Stability Analysis
r = roots([1;-1*Wc]);
fprintf('The roots of characteristic equation are:\n');
disp(r);
fprintf('The magnitudes of roots are: \n');
disp(abs(r));

% Goodness-of-fit Metrics
mse = mean(e(wlen+1:N).^2); % Mean-squared error
REV = mse/var(xfn(wlen+1:N)); % Relative Error Variance
% Normalized Residual
norm_res = sqrt(sum(e(wlen+1:N,1).^2))/sqrt(sum(xfn(wlen+1:N,1).^2));
fprintf('Relative Error Variance is %d\n',REV);
fprintf('Normalized Residual is %d\n',norm_res);

%% Display Results
figure('color',[1 1 1]);
subplot(2,1,1);
plot(xfn(wlen+1:N),'r');
hold on;
plot(x_est(wlen+1:N),'--b','Linewidth',1);
ylim([-1.2,1.2]);
xlabel('Time (ms)','Fontsize',12);
ylabel('Normalized Amplitude','Fontsize',12);
title(' A Input vs Predicted LFPs','Fontweight','b','Fontsize',12);
legend('Input LFP','Predicted LFP');

figure('color',[1 1 1]);
subplot(2,2,1);
plot(x);
ax1 = gca;
ax1.FontSize = 10;
ax1.FontWeight = 'b';
ax1.YTick = [-1e-3:0.5e-3:1e-3];
ax1.YTickLabel = {[-1:0.5:1]};
ylim([-1e-3,1e-3]);
title(' A Raw LFP','Fontsize',12);
xlabel('Time (ms)','Fontsize',12);
ylabel('Amplitude (mV)','Fontsize',12);

subplot(2,2,2);
plot(f_fft,x_psd);
ax1 = gca;
ax1.FontSize = 10;
ax1.FontWeight = 'b';
ax1.XTick = [3,5,7,9,11,13,15,20];
xlim([3,15]);
title(' B Fourier Spectrum of Raw LFP','Fontsize',12);
xlabel('Frequency (Hz)','Fontsize',12);
ylabel('PSD (V^2/Hz)','Fontsize',12);

subplot(2,2,3);
plot(wt(1,:),'r');
hold on;
plot(wt(2,:),'g');
hold on;
plot(wt(3,:),'k');
hold on;
plot(wt(4,:),'b');
hold on;
plot(wt(5,:),'c');
hold on;
plot(wt(6,:),'m');
legend('w(1)','w(2)','w(3)','w(4)','w(5)','w(6)');
ax1 = gca;
ax1.FontSize = 10;
ax1.FontWeight = 'b';
title(' C Trajectory of AR Model Coefficients','Fontsize',12);
xlabel('Time (ms)','Fontsize',12);
ylabel('Amplitude','Fontsize',12);

subplot(2,2,4);
% plot(f,MinMaxNorm(h_mag,0,1));
plot(f,h_mag);
ax1 = gca;
ax1.FontSize = 10;
ax1.FontWeight = 'b';
ax1.XTick = [3,5,7,9,11,13,15,20];
xlim([3,15]);
title(' D Frequency Response of Offline-Learned AR Coefficients','Fontsize',12);
xlabel('Frequency (Hz)','Fontsize',12);
ylabel('Magnitude','Fontsize',12);
