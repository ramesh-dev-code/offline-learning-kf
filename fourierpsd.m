function [f,psd] = fourierpsd(X,fs)
L = length(X); % Length of input signal
N = 2^nextpow2(L);
f = 0:fs/N:fs/2; % Frequency points in PSD plot
Xf = fft(X,N);
psd = (abs(Xf(1:N/2+1)).^2)/N; % power spectral distribution
psd(2:end-1) = 2*psd(2:end-1); % single-sided psd 
% figure('color',[1 1 1]);
% subplot(2,2,1);
% plot(X,'Linewidth',2);
% ylim([-1e-6,1e-6]);
% title('  (A) Input Signal','Fontsize',12,'Fontweight','b');
% xlabel('Time (ms)','Fontsize',12,'Fontweight','b');
% ylabel('Voltage (\muV)','Fontsize',12,'Fontweight','b');
% subplot(2,2,2);
% plot(f,psd,'Linewidth',2);
% xlim([0,50]);
% % ylim([0,10e-17]);
% title('  (B) Fourier Power Spectrum','Fontsize',14,'Fontweight','b');
% xlabel('Frequency (Hz)','Fontsize',12,'Fontweight','b');
% ylabel('PSD','Fontsize',12,'Fontweight','b');
end