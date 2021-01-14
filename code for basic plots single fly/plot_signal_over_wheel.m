function [out] = plot_signal_over_wheel(recording, ROI_num)


signal_store = [];
signal_time_store = [];
filtered_wheel_store = [];
sucrose_store = [];

if(ROI_num == -1)
    recording.tseries = 1;
end
for i = 1:1:length(recording.tseries)
    if(ROI_num == -1)
     signal_time = recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
        
         signal    = recording.abf.CH1_patch(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000))-13;
         signal    = recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000));
    else
    signal      = recording.tseries(i).df_over_f(ROI_num,:);
    signal_time = recording.tseries(i).Time_s;
    end
    interp_filt_wheel = interp1(recording.movie1.time_stamps, recording.movie1.filtered_wheel, signal_time,'previous');
    interp_sucrose = interp1(recording.movie1.time_stamps, recording.movie1.sucrose, signal_time,'previous');
    
    signal_store = [signal_store; signal'];
    signal_time_store = [signal_time_store; signal_time];
    filtered_wheel_store = [filtered_wheel_store; interp_filt_wheel];
    sucrose_store = [sucrose_store, interp_sucrose'];
end

figure; hold on;
set(gca,'TickDir','out');

%subplot(1,2,1); hold on;
title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
[a b] = find(sucrose_store' ==200);
scatter(filtered_wheel_store(a), signal_store(a), 3.*ones(1,length(signal_store(a))), [33 113 181]./255, 'filled');
mean_sucrose200 = nanmean(signal_store(a));
sem_sucrose200 = nanstd(signal_store(a))./sqrt(length(a));


[a b] = find(sucrose_store' ==500);
scatter(filtered_wheel_store(a), signal_store(a), 3.*ones(1,length(signal_store(a))), [8,48,107]./255 , 'filled');
mean_sucrose500 = nanmean(signal_store(a));
sem_sucrose500 = nanstd(signal_store(a))./sqrt(length(a));

[a b] = find(sucrose_store' == 0);
hold on; scatter(filtered_wheel_store(a), signal_store(a), 3.*ones(1,length(signal_store(a))), [107,174,214]./255, 'filled');
mean_sucrose0 = nanmean(signal_store(a));
sem_sucrose0 = nanstd(signal_store(a))./sqrt(length(a));
cnt = 1;
m = [];
medin = [];

tmp_filtered_wheel_store = [filtered_wheel_store; filtered_wheel_store+2*pi; filtered_wheel_store+4*pi; filtered_wheel_store+6*pi];
tmp_signal_store = [signal_store; signal_store; signal_store; signal_store];

for iterat = 0:pi/50:(6*pi)
    [a b] = find(tmp_filtered_wheel_store >= iterat & tmp_filtered_wheel_store < (iterat+pi/4));
    m(cnt) = nanmean(tmp_signal_store(a));
    medin(cnt) = nanmedian(tmp_signal_store(a));
    
    cnt = cnt+1;
end

[a b] = sort(mod(pi/8+(0:pi/50:6*pi),2*pi),'ascend');
plot(a, m(b),'k','linewidth',2);
plot(a, medin(b),'Color',[.4 .4 .4],'linewidth',2);
xlim([0,2*pi]);

figure; hold on;
title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');

set(gca,'TickDir','out');

%subplot(1,2,2);  hold on;
errorbar(0,mean_sucrose500,sem_sucrose500,'o','LineWidth',2,'Color',[8,48,107]./255);
text(0,mean_sucrose500+.3,num2str(mean_sucrose500));

errorbar(1,mean_sucrose200,sem_sucrose200,'o','LineWidth',2,'Color',[33 113 181]./255);
text(1,mean_sucrose200+.3,num2str(mean_sucrose200));

errorbar(2,mean_sucrose0,sem_sucrose0,'o','LineWidth',2,'Color',[107,174,214]./255);
text(2,mean_sucrose0+.3,num2str(mean_sucrose0));



xticks([0 1 2])
xticklabels({'500', '200', '0'})


xlim([-2,4]);
ylim([-.5, 3]);



out=1;

end

