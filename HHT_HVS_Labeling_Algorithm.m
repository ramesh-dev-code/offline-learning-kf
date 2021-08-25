close all;
clear all;
clc;
dataset_type = 'testing';
dataset_name = 'R1_M1D';
data = importdata(char(strcat('GroundTruth\',dataset_type,'\',dataset_name,'.csv')));
x = data(:,1);
N = length(x);
%%
win_t = 200;
xn = varnorm1(x,0,1);
t1 = 57001;
t2 = t1 + 999;
fs = 1000;
fmax = 50; 
f_res = 0.5; % Unit: Hz
fl = 11;
fh = 28;
[imf,residual,info] = emd(x);
num_imf = size(imf,2);
% hs - time-frequency distribution of the Hilbert transform envelope amplitude 
[hs,f,t,imf_ins_freq,imf_ins_energy] = hht(imf,fs,'FrequencyLimits',[0,fmax],'FrequencyResolution',f_res);
m_hvs_ie = mean(full(hs(fl:fh,:)),1); % Inst Mean HVS Energy in 5~13 Hz band 
m_hvs_ie = filter(gausswin(win_t),1,m_hvs_ie); % Windowed Mean HVS Energy
%%
beta = 10; % Scaling Factor
% Uncomment the following line to compute the detection threshold from the
% validation dataset
% tr2 = mean(m_hvs_ie(t1:t2)) + beta*std(m_hvs_ie(t1:t2));
tr2 = 2.1e-8;   
clear hvs3;
hvs3 = m_hvs_ie>tr2;  
hvs4 = hvs3;
CL = 200;
t_st = 1;
if CL > 0
    for i = 1:N-CL
        % Remove false positives
        if hvs3(i)==1 && ~all(hvs3(i+1:i+CL))     
            hvs4(i) = 0;
        end    
    end
end

% Detection response based on m_hvs_ie
fig = figure('Color',[1 1 1]);
subplot(2,1,1);
plot(t_st:N,xn(t_st:N),'r');
hold on;
plot(t_st:N,hvs4(t_st:N),'k','Linewidth',1.5);
xlim([t_st,N]);
ylim([-0.2,1.2]);
title('A Detection Response','FontSize',10,'FontWeight','b');
xlabel('Time (ms)','FontSize',10,'FontWeight','b');
ylabel('Normalized Amplitude','FontSize',10,'FontWeight','b');
datacursormode on
dcm_obj = datacursormode(fig);
% Set update function
set(dcm_obj,'UpdateFcn',@myupdatefcn)
pause 
% Export cursor to workspace
info_struct = getCursorInfo(dcm_obj);
if isfield(info_struct, 'Position')
  disp('Clicked positioin is')
  disp(info_struct.Position)
end

subplot(2,1,2);
plot(t_st:N,m_hvs_ie(t_st:N),'k','Linewidth',1.5);
title('B Instantaneous Mean Energy in HVS Band','FontSize',10,'FontWeight','b');
xlabel('Time (ms)','FontSize',10,'FontWeight','b');
ylabel('Mean HVS Energy (V^2)','FontSize',10,'FontWeight','b');
xlim([t_st,N]);
ylim([0,1e-6]);

%% Uncomment this segment after manually correcting the false positives
% clear hvs;
% clear onset;
% clear hvs_end;
% fname = char(strcat('GroundTruth\',dataset_type,'\',dataset_name,'.csv'));
% hvs = diff([0,hvs4]);
% onset = [find(hvs>0)]';
% hvs_end = [find(hvs<0)-1]';
% csvwrite(fname,[x,m_hvs_ie',hvs4']);
