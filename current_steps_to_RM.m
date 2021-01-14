%% This code takes care of Ra correction and -13 mV correction
%% Also, takes care of bias (usually -2 pA) in the current recorded in abf. Even at I=0 the abf shows a bit of current as being injected. This is corrected for.



% mean_for_individual_flies = [];
% mean_for_individual_flies.Vave = [];
% mean_for_individual_flies.Iave = [];
% mean_for_individual_flies.Spikeave = [];
% cnt = 1;


%% VT kir
%% abfs from an individual cell (repeats)
% 
abfs = {'2020_10_08_0003.abf','2020_10_08_0004.abf'};
Ra = 17.3e6;

 abfs = {'2020_10_08_0006.abf','2020_10_08_0006.abf'};
 Ra = 20e6;
%  
abfs = {'2020_10_08_0009.abf','2020_10_08_0010.abf'};
Ra = 13.7e6;
%  
abfs = {'2020_10_08_0013.abf','2020_10_08_0013.abf'};
Ra = 12.5e6;
 
% same cell as above, different current steps
% abfs = {'2020_10_08_0014.abf','2020_10_08_0014.abf'};
% Ra = 12.5e6;
 
% same cell as above, different current steps
% abfs = {'2020_10_08_0015.abf','2020_10_08_0017.abf'};
% Ra = 12.5e6;
 
abfs = {'2020_10_08_0019.abf','2020_10_08_0020.abf','2020_10_08_0021.abf','2020_10_08_0022.abf','2020_10_08_0023.abf'};
Ra = 12.9e6;
 
%% VT mutant
%% abfs from an individual cell (repeats)

abfs = {'2020_10_09_0002.abf','2020_10_09_0003.abf','2020_10_09_0004.abf'};
Ra = 22.7e6;

abfs = {'2020_10_13_0001.abf','2020_10_13_0002.abf','2020_10_13_0003.abf'};
Ra = 14.2e6;

abfs = {'2020_10_13_0005.abf','2020_10_13_0006.abf','2020_10_13_0007.abf'};
Ra = 12.1e6;
% 
% 
abfs = {'2020_10_13_0010.abf','2020_10_13_0011.abf','2020_10_13_0012.abf'};
Ra = 2.1e6;

%% SS kir and mutant

%% 2xegfp
% this recording was before I had current in abf
% abfs = {'2020_09_09_0002.abf','2020_09_09_0003.abf','2020_09_09_0004.abf'};
% Ra = 24.5e6;
 


% .2 to 1.05 is stable for the 5 sec sweeps
interval = 2000:1:10500;
interval_baseline = 30000:1:50000;

% interval = 1000:1:5000;
% interval_baseline = 6500:1:10000;

% Rm to average (here it is set from -15 to -40 pA)
Rm_to_ave = [-15,-35];






V_ave = [];
I_ave = [];
R_ave = [];

V_median = [];
I_median = [];
R_median = [];
spikes_per_sec = [];

V_all = [];
I_all = [];
R_all = [];

V_baseline = [];
I_baseline = [];

   f1 = figure; hold on;

for j = 1:1:length(abfs)
    
    q = cell2str(abfs(j));
    q = q(4:(end-3));
    
    [d,si,h]=abfload(q);
    
% only required for the 09/09/2020 recording that dosent have current in
% the abf
% d(interval_baseline,2,:) = 0;
% steps = -50:5:45;
% for i = 1:1:20
%     d(interval,2,i) = steps(i);
% end

    
    V_baseline(j) = mean(mean(d(interval_baseline,1,:)))-13;
    I_baseline(j) = mean(mean(d(interval_baseline,2,:)));

    
    for i = 1:1:length(d(1,1,:))
        I_ave(j,i) = mean(d(interval,2,i))-I_baseline(j);
        V_ave(j,i) = mean(d(interval,1,i))-13-(Ra*I_ave(j,i)*1e-12)*1000;
        R_ave(j,i) = V_ave(i)/I_ave(i);
        
        I_all = [I_all; d(interval,2,i)-I_baseline(j)];
        V_all = [V_all; d(interval,1,i)-13-(Ra.*d(interval,2,i).*1e-12).*1000];
        
        
        I_median(j,i) = median(d(interval,2,i))-I_baseline(j);
        V_median(j,i) = median(d(interval,1,i))-13-(Ra*I_median(j,i)*1e-12)*1000;
        R_median(j,i) = V_median(i)/I_median(i);
        
        [spikes_out]=calculate_spikes(d(interval,1,i),100,1.4);
        spikes_out_per_second = sum(spikes_out)/ ((interval(end)-interval(1))./10000);
        spikes_per_sec(j,i) = spikes_out_per_second;

        d(:,3,i) = d(:,1,i) -13-(Ra*d(:,2,i)*1e-12)*1000;
    end
    
    
    
 figure(f1); hold on;
 set(gca,'ColorOrderIndex',1)
for i = 1:1:length(d(1,3,:))
plot(d(:,3,i));
end

ylabel('mean corrected Vm during interval (mV)');
xlabel('samples (10,000/sec)');
set(gca,'TickDir','out')
set(gca,'xlim',[0,20000]);
set(gca,'ylim',[-120,10]);


end
R_all = V_all./I_all;

mean_for_individual_flies(cnt).Vave = mean(V_ave);
mean_for_individual_flies(cnt).Iave = mean(I_ave);
mean_for_individual_flies(cnt).Spikeave = mean(spikes_per_sec);
cnt = cnt+1;

% figure; hold on;title('');
% ylabel('spikes per second');
% xlabel('mean corrected Vm  during interval (mV)');
% set(gca,'TickDir','out')
% %line([-60,60],[mean(V_baseline),mean(V_baseline)],'color','k');
% %scatter(V_ave,spikes_per_sec);
% [size_d1 size_d2] = size(I_ave);
% plot(V_ave,spikes_per_sec,'Markersize',3,'Marker','o','linestyle','none');
% errorbar(mean(V_ave),mean(spikes_per_sec),std(spikes_per_sec)./sqrt(size_d1),std(spikes_per_sec)./sqrt(size_d1),std(V_ave)./sqrt(size_d1),std(V_ave)./sqrt(size_d1),'k');
% set(gca,'ylim',[0,25]);
% set(gca,'xlim',[-100,10]);
% set(gca,'TickDir','out')
% 
% figure; hold on;title('');
% ylabel('spikes per second');
% xlabel('mean current during interval');
% set(gca,'TickDir','out')
% %line([-60,60],[mean(V_baseline),mean(V_baseline)],'color','k');
% %scatter(I_ave,spikes_per_sec);
% [size_d1 size_d2] = size(I_ave);
% plot(I_ave,spikes_per_sec,'Markersize',3,'Marker','o','linestyle','none');
% errorbar(mean(I_ave),mean(spikes_per_sec),std(spikes_per_sec)./sqrt(size_d1),std(spikes_per_sec)./sqrt(size_d1),std(I_ave)./sqrt(size_d1),std(I_ave)./sqrt(size_d1),'k');
% set(gca,'ylim',[0,25]);
% set(gca,'xlim',[-50,50]);
% set(gca,'TickDir','out')








% figure; hold on;title('slope is Rm (mean during pulse)');
% ylabel('mean corrected Vm during interval (mV)');
% xlabel('mean current during interval (pA)');
% set(gca,'TickDir','out')
% line([-60,60],[mean(V_baseline),mean(V_baseline)],'color','k');
% %scatter(V_ave,I_ave);
% plot(I_ave,V_ave,'Marker','o');
% set(gca,'ylim',[-150,0]);
% set(gca,'xlim',[-60,60]);

% figure; hold on; title('slope is Rm (median during pulse)');
% ylabel('median corrected Vm during interval (mV)');
% xlabel('median current during interval (pA)');
% line([-60,60],[mean(V_baseline),mean(V_baseline)],'color','k');
% set(gca,'TickDir','out')
% %scatter(V_median,I_median);
% plot(I_median,V_median,'Marker','o');
% set(gca,'ylim',[-150,0]);
% set(gca,'xlim',[-60,60]);

% figure; hold on; title('Rm (mean during pulse)');
% xlabel('mean corrected Vm during interval (mV)');
% ylabel('Rm (MOhm)');
% line([mean(V_baseline),mean(V_baseline)],[-100,5000],'color','k','Linestyle','--');
% set(gca,'TickDir','out')
% %scatter(V_median,I_median);
% plot(V_ave,1e-6*(V_ave./1000-mean(V_baseline)./1000)./(I_ave.*1e-12),'Marker','o','Linestyle','none');
% errorbar(mean(V_ave),1e-6*mean((V_ave./1000-mean(V_baseline)./1000)./(I_ave.*1e-12)),1e-6*std((V_ave./1000-mean(V_baseline)./1000)./(I_ave.*1e-12)),'k','Linestyle','none');
% tmp1 = mean((V_ave./1000-mean(V_baseline)./1000)./(I_ave.*1e-12));
% tmp2 = mean(V_ave);
% tmp3 = mean((I_ave));
% 
% for ii = 1:1:length(mean(V_ave))
%     h = text(tmp2(ii), 2800, num2str(round(tmp1(ii)./1e6)),'Fontsize',7);
%     set(h,'Rotation',90);
% end
% set(gca,'ylim',[0,3000]);
% set(gca,'xlim',[-150,0]);
% 
% [a b] = find(tmp3 <= Rm_to_ave(1) & tmp3 >= Rm_to_ave(2));
%     h = text(-100, 200, num2str(mean(round(tmp1(b)./1e6))),'Fontsize',12);


% 
% figure; hold on; title('slope is Rm (individual samples)');
% ylabel('mean corrected Vm during interval (mV)');
% xlabel('mean current during interval (pA)');
% line([-60,60],[mean(V_baseline),mean(V_baseline)],'color','k');
% set(gca,'TickDir','out')
% scatter(I_all,V_all,3,'filled','k');
% set(gca,'ylim',[-150,0]);
% set(gca,'xlim',[-60,60]);

