function [out] = plot_signal_at_transitions(recording, ROI_num, filter_out_chrimson_data)

%% plot signal after transitions

% this is all transitions and never plots a second transitions
% data is stopped when there is another transition in the positive or
% negative portions of this plot
% set ROI num to -1 if you want tu use CH1 patch
% if patch is 1 then use CH1
% if patch is 0 then use CH1 but remove pulses (and 100 seconds after)

  if(ROI_num ==-1)
        patch = 1;
        recording.tseries = [];
    else
        patch = 0;
  end

cntr = 0;
cell_store_trans = {};
cell_store_trans_timebase = {};
cell_store_filtered_wheel = {};
cell_store_sucrose = {};
store_trans = [];

out_signal = [];
out_sucrose = [];
out_filtered_wheel = [];
out_timebase = [];

for ts = 1:1:length(recording.tseries)
    ds = diff(recording.movie1.sucrose);
    ds = [0,ds];
    [a b] = find(ds);
    b = [1,b, length(recording.movie1.sucrose)];
    
    for i =2:1:(length(b)-1)
        
        startt =recording.movie1.time_stamps(b(i-1));
        stopt = recording.movie1.time_stamps(b(i+1)-1);
        
        if(patch ==0)
            [time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num,:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
        end
        
        if(patch ==1)
            [time_base_to_return, data_to_average_interp] = average_around_event(recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
        end
        
        if(patch ==0 && filter_out_chrimson_data == 1)
            
            tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(ts).Time_s,'previous')';
            [time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num,:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2, recording.abf.Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
            data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
            data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
        end
        
        if(patch ==1  && filter_out_chrimson_data == 1)
            [time_base_to_return, data_to_average_interp] = average_around_event(recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
            data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
            data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
            
        end
        
        
        
        
        
        
        
        
        
        [time_base_to_return1, data_to_average_interp1] = average_around_event(recording.movie1.sucrose, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
        [time_base_to_return2, data_to_average_interp2] = average_around_event(recording.movie1.filtered_wheel, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
        
        % time between transition (i-1) and i+1 has to be greater
        % than 4 seconds and has to be in transition zone (removes
        % errant)
        
        if( (b(i+1)-b(i-1)) > 100)
            
            cntr = cntr+1;
            
            cell_store_trans{cntr} =  data_to_average_interp;
            cell_store_trans_timebase{cntr} = time_base_to_return;
            cell_store_filtered_wheel{cntr} = data_to_average_interp2;
            cell_store_sucrose{cntr} = data_to_average_interp1;
            store_trans(cntr,:) = [b(i-1), b(i), b(i+1)-1];
        end
    end
end

for i = 1:1:cntr
    out_signal(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans{i}, -240:.1:900,'previous');
    out_sucrose(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_sucrose{i}, -240:.1:900,'previous');
    out_filtered_wheel(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_filtered_wheel{i}, -240:.1:900,'previous');
    out_timebase(i,:) = interp1(cell_store_trans_timebase{i}, cell_store_trans_timebase{i}, -240:.1:900,'previous');
end

out_signal_copy = out_signal;
out_signal(~any(~isnan(out_signal), 2),:)=[];
out_sucrose(~any(~isnan(out_signal_copy), 2),:)=[];
out_filtered_wheel(~any(~isnan(out_signal_copy), 2),:)=[];
out_timebase(~any(~isnan(out_signal_copy), 2),:)=[];
store_trans(~any(~isnan(out_signal_copy), 2),:)=[];

% plotting after calculations are made
% no filtering of data here (except transition time)
% first just plot where each of the transitions are

%%%%%%



[rr,c] = size(out_signal);
s200_to_p = [];
p_to_s200 = [];

s200_to_s500 = [];
s500_to_s200 = [];

s500_to_p = [];
p_to_s500 = [];

f_200_to_0 = figure;
f_0_to_200 = figure;
f_200_to_500 = figure;
f_500_to_200 = figure;
f_500_to_0 = figure;
f_0_to_500 = figure;


for r = 1:1:rr
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
        figure(f_0_to_200); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) == 0) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 1 .75]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) == 200) = NaN;
        %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        p_to_s200 = [p_to_s200, r];
    end
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
        figure(f_200_to_0); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 200) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 1 .75]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 0) = NaN;
        %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        s200_to_p = [s200_to_p, r];
        
    end
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
        figure(f_0_to_500); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 500) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 .75 1]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 0) = NaN;
        %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        p_to_s500 = [p_to_s500, r];
    end
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
        figure(f_500_to_0); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 500) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 .75 1]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 0) = NaN;
        %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        s500_to_p = [s500_to_p, r];
        
    end
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
        figure(f_200_to_500); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 500) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 .75 1]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 200) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 1  .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        s200_to_s500 = [s200_to_s500, r];
    end
    
    if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
        figure(f_500_to_200); hold on;
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 500) = NaN;
        %plot(-240:.1:900,qq,'color',[.75 .75 1]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        qq = out_signal(r,:);
        qq(out_sucrose(r,:) ~= 200) = NaN;
        %    plot(-240:.1:900,qq,'color',[.75 1 .75 ]);
        plot(-240:.1:900,qq,'color',[220 220 220]./255);
        
        s500_to_s200 = [s500_to_s200, r];
        
    end
    
end



figure(f_200_to_0); hold on; title(['200 mM sucrose to plain transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(s200_to_p))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_0_to_200); hold on; title(['plain to 200 mM sucrose transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(p_to_s200))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_500_to_0); hold on; title(['500 mM sucrose to plain transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(s500_to_p))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_0_to_500); hold on; title(['plain to 500 mM sucrose transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(p_to_s500))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_200_to_500); hold on; title(['200 mM sucrose to 500 mM sucrose transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(s200_to_s500))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_500_to_200); hold on; title(['500 mM sucrose to 200 mM sucrose transitions' ' ' 'ROI ' num2str(ROI_num) ' total fragments ' num2str(length(s500_to_s200))]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

analysis = figure;

%%%%
figure(f_200_to_0); hold on;
qq = out_signal(s200_to_p,:);
qq(out_sucrose(s200_to_p,:) ~= 200) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','g');
plot(-240:.1:900,nanmean(qq),'color',[33 113 181]./255,'linewidth',2);

% 
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% 
% mean_sucrose200 = nanmean(qq);
% sem_sucrose200 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(s200_to_p,:);
qq(out_sucrose(s200_to_p,:) ~= 0) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','r');
plot(-240:.1:900,nanmean(qq),'color',[107 174 214]./255,'linewidth',2);
% 
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% 
% mean_sucrose0 = nanmean(qq);
% sem_sucrose0 = nanstd(qq)./sum(~isnan(qq));

xlim([-45,45]);

% figure(analysis); hold on;
% 
% errorbar(0,mean_sucrose200,sem_sucrose200,'o','LineWidth',2,'Color',[33 113 181]./255);
% text(0,mean_sucrose200+.3,num2str(mean_sucrose200));
% 
% errorbar(1,mean_sucrose0,sem_sucrose0,'o','LineWidth',2,'Color',[107,174,214]./255);
% text(1,mean_sucrose0+.3,num2str(mean_sucrose0));

% 
% text(.5,2,num2str(100*( mean_sucrose200-mean_sucrose0) / mean_sucrose0));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(f_0_to_200); hold on;
qq = out_signal(p_to_s200,:);
qq(out_sucrose(p_to_s200,:) ~= 200) = NaN;
% plot(-240:.1:900,nanmean(qq),'color','g');
plot(-240:.1:900,nanmean(qq),'color',[33 113 181]./255,'linewidth',2);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose200 = nanmean(qq);
% sem_sucrose200 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(p_to_s200,:);
qq(out_sucrose(p_to_s200,:) ~= 0) = NaN;
% plot(-240:.1:900,nanmean(qq),'color','r');
plot(-240:.1:900,nanmean(qq),'color',[107 174 214]./255,'linewidth',2);


xlim([-45,45]);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose0 = nanmean(qq);
% sem_sucrose0 = nanstd(qq)./sum(~isnan(qq));


figure(analysis); hold on;
set(gca,'TickDir','out');

% errorbar(4,mean_sucrose200,sem_sucrose200,'o','LineWidth',2,'Color',[33 113 181]./255);
% text(4,mean_sucrose200+.3,num2str(mean_sucrose200));
% 
% errorbar(3,mean_sucrose0,sem_sucrose0,'o','LineWidth',2,'Color',[107,174,214]./255);
% text(3,mean_sucrose0+.3,num2str(mean_sucrose0));
% 
% text(3.5,2,num2str(100*(mean_sucrose0-mean_sucrose200)/mean_sucrose200));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
figure(f_500_to_0); hold on;
qq = out_signal(s500_to_p,:);
qq(out_sucrose(s500_to_p,:) ~= 500) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','b');
plot(-240:.1:900,nanmean(qq),'color',[8 48 107]./255,'linewidth',2);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose500 = nanmean(qq);
% sem_sucrose500 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(s500_to_p,:);
qq(out_sucrose(s500_to_p,:) ~= 0) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','r');
plot(-240:.1:900,nanmean(qq),'color',[107 174 214]./255,'linewidth',2);


xlim([-45,45]);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose0 = nanmean(qq);
% sem_sucrose0 = nanstd(qq)./sum(~isnan(qq));


% figure(analysis); hold on;
% 
% errorbar(6,mean_sucrose500,sem_sucrose500,'o','LineWidth',2,'Color',[8 48 107]./255);
% text(6,mean_sucrose500+.3,num2str(mean_sucrose500));
% 
% errorbar(7,mean_sucrose0,sem_sucrose0,'o','LineWidth',2,'Color',[107,174,214]./255);
% text(7,mean_sucrose0+.3,num2str(mean_sucrose0));
% 
% text(6.5,2,num2str(100*(mean_sucrose500-mean_sucrose0)/mean_sucrose0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(f_0_to_500); hold on;
qq = out_signal(p_to_s500,:);
qq(out_sucrose(p_to_s500,:) ~= 500) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','b');
plot(-240:.1:900,nanmean(qq),'color',[8 48 107]./255,'linewidth',2);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose500 = nanmean(qq);
% sem_sucrose500 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(p_to_s500,:);
qq(out_sucrose(p_to_s500,:) ~= 0) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','r');
plot(-240:.1:900,nanmean(qq),'color',[107 174 214]./255,'linewidth',2);


xlim([-45,45]);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose0 = nanmean(qq);
% sem_sucrose0 = nanstd(qq)./sum(~isnan(qq));


% figure(analysis); hold on;
% 
% errorbar(10,mean_sucrose500,sem_sucrose500,'o','LineWidth',2,'Color',[8 48 107]./255);
% text(10,mean_sucrose500+.3,num2str(mean_sucrose500));
% 
% errorbar(9,mean_sucrose0,sem_sucrose0,'o','LineWidth',2,'Color',[107,174,214]./255);
% text(9,mean_sucrose0+.3,num2str(mean_sucrose0));
% 
% 
% text(9.5,2,num2str(100*(mean_sucrose0-mean_sucrose500)/mean_sucrose500));


%%%%
figure(f_500_to_200); hold on;
qq = out_signal(s500_to_s200,:);
qq(out_sucrose(s500_to_s200,:) ~= 500) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','b');
plot(-240:.1:900,nanmean(qq),'color',[8 48 107]./255,'linewidth',2);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose500 = nanmean(qq);
% sem_sucrose500 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(s500_to_s200,:);
qq(out_sucrose(s500_to_s200,:) ~= 200) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','g');
plot(-240:.1:900,nanmean(qq),'color',[33 113 181]./255,'linewidth',2);


xlim([-45,45]);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose200 = nanmean(qq);
% sem_sucrose200 = nanstd(qq)./sum(~isnan(qq));

% 
% figure(analysis); hold on;
% 
% errorbar(12,mean_sucrose500,sem_sucrose500,'o','LineWidth',2,'Color',[8 48 107]./255);
% text(12,mean_sucrose500+.3,num2str(mean_sucrose500));
% 
% errorbar(13,mean_sucrose200,sem_sucrose200,'o','LineWidth',2,'Color',[33 113 181]./255);
% text(13,mean_sucrose200+.3,num2str(mean_sucrose200));
% 
% text(12.5,2,num2str(100*(mean_sucrose500-mean_sucrose200)/mean_sucrose200));
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(f_200_to_500); hold on;
qq = out_signal(s200_to_s500,:);
qq(out_sucrose(s200_to_s500,:) ~= 500) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','b');
plot(-240:.1:900,nanmean(qq),'color',[8 48 107]./255,'linewidth',2);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose500 = nanmean(qq);
% sem_sucrose500 = nanstd(qq)./sum(~isnan(qq));

qq = out_signal(s200_to_s500,:);
qq(out_sucrose(s200_to_s500,:) ~= 200) = NaN;
%plot(-240:.1:900,nanmean(qq),'color','g');
plot(-240:.1:900,nanmean(qq),'color',[33 113 181]./255,'linewidth',2);


xlim([-45,45]);
% qq = qq(:,1950:(1950+900));
% qq = qq(:);
% mean_sucrose200 = nanmean(qq);
% sem_sucrose200 = nanstd(qq)./sum(~isnan(qq));


% figure(analysis); hold on;
% 
% errorbar(16,mean_sucrose500,sem_sucrose500,'o','LineWidth',2,'Color',[8 48 107]./255);
% text(16,mean_sucrose500+.3,num2str(mean_sucrose500));
% 
% errorbar(15,mean_sucrose200,sem_sucrose200,'o','LineWidth',2,'Color',[33 113 181]./255);
% text(15,mean_sucrose200+.3,num2str(mean_sucrose200));
% 
% text(15.5,2,num2str(100*(mean_sucrose200-mean_sucrose500)/mean_sucrose500));

% xlim([-2,17]);
% ylim([-.5, 3]);
% 
% xticks([.5 3.5 6.5 9.5 12.5 15.5])
% xticklabels({'200->0','0->200','500->0','0->500','500->200','200->500'})
% 
% snapnow;







figure(f_200_to_0); snapnow;

figure(f_0_to_200); snapnow;

figure(f_500_to_0); snapnow;

figure(f_0_to_500); snapnow;

figure(f_200_to_500); snapnow;

figure(f_500_to_200); snapnow;


figure(f_200_to_0); xlim([-120,120]);   snapnow;


figure(f_0_to_200); xlim([-120,120]);  snapnow;
figure(f_500_to_0); xlim([-120,120]);   snapnow;

figure(f_0_to_500); xlim([-120,120]);  snapnow;

figure(f_200_to_500); xlim([-120,120]);  snapnow;

figure(f_500_to_200); xlim([-120,120]);   snapnow;

out=1;

end

