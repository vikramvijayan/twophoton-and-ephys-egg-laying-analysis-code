%egg_structure = egg_structure_median;



patch = 0;
time_base_egg = -240:.1:240;

% plotting egg triggered averages
data_around_egg_vel = [egg_structure(1:end).data_around_egg_vel];
data_around_egg_vel_noave = [egg_structure(1:end).data_around_egg_vel_noave];

for i =1:1:length(egg_structure)
    egg_structure(i).data_around_egg_sub = egg_structure(i).data_around_egg_sub';
end
data_around_egg_sub =   [egg_structure(1:end).data_around_egg_sub];

data_around_egg_body = [egg_structure(1:end).data_around_egg_body];
data_around_egg_body_x = [egg_structure(1:end).data_around_egg_body_x];
data_around_egg_body_y = [egg_structure(1:end).data_around_egg_body_y];

data_around_egg_prob = [egg_structure(1:end).data_around_egg_prob];
data_around_egg_prob_x = [egg_structure(1:end).data_around_egg_prob_x];
data_around_egg_prob_y = [egg_structure(1:end).data_around_egg_prob_y];

data_around_egg_path = [egg_structure(1:end).data_around_egg_path];
data_around_egg_body_only = [egg_structure(1:end).data_around_egg_body_only];
data_around_egg_body_only_x = [egg_structure(1:end).data_around_egg_body_only_x];
data_around_egg_body_only_y = [egg_structure(1:end).data_around_egg_body_only_y];

data_around_egg_angle = [egg_structure(1:end).data_around_egg_angle];
data_around_egg_angle_only = [egg_structure(1:end).data_around_egg_angle_only];
data_around_egg = [egg_structure(1:end).data_around_egg];


% remove double plotting of behavior
if(~isempty(data_around_egg))
    [a b] = find(~isnan(data_around_egg(2400,:)));
    
    data_around_egg_vel = data_around_egg_vel(:,b);
    
    data_around_egg_vel_noave = data_around_egg_vel_noave(:,b);
    
    data_around_egg_body = data_around_egg_body(:,b);
    data_around_egg_body_x = data_around_egg_body_x(:,b);
    data_around_egg_body_y = data_around_egg_body_y(:,b);
    
    data_around_egg_prob = data_around_egg_prob(:,b);
    data_around_egg_prob_x = data_around_egg_prob_x(:,b);
    data_around_egg_prob_y = data_around_egg_prob_y(:,b);
    
    data_around_egg_path = data_around_egg_path(:,b);
    data_around_egg_body_only = data_around_egg_body_only(:,b);
    data_around_egg_body_only_x = data_around_egg_body_only_x(:,b);
    data_around_egg_body_only_y =data_around_egg_body_only_y(:,b);
    
    data_around_egg_angle = data_around_egg_angle(:,b);
    data_around_egg_angle_only = data_around_egg_angle_only(:,b);
    data_around_egg = data_around_egg(:,b);
    data_around_egg_sub = data_around_egg_sub(:,b);
end

[pp,ppp] = size(data_around_egg_body);
time_base_pulse = [-240:.1:240];

% flter out large jumps in velocity and proboscis
data_around_egg_vel_noave_filt = data_around_egg_vel_noave;

[atemp_diff_vel btemp_diff_vel] = size(data_around_egg_vel_noave);
for vel_col = 1:1:btemp_diff_vel
    for vel_row = 2:1:atemp_diff_vel
        if(vel_row == 2)
            data_around_egg_vel_noave_filt(vel_row-1,vel_col) = nanmedian(data_around_egg_vel_noave_filt(:,vel_col));
        end
        if(abs(data_around_egg_vel_noave_filt(vel_row,vel_col)-data_around_egg_vel_noave_filt(vel_row-1,vel_col)) > 10)
            data_around_egg_vel_noave_filt(vel_row,vel_col) = data_around_egg_vel_noave_filt(vel_row-1,vel_col);
        end
    end
end




data_around_egg_vel_ave_filt = data_around_egg_vel;

[atemp_diff_vel btemp_diff_vel] = size(data_around_egg_vel);
for vel_col = 1:1:btemp_diff_vel
    for vel_row = 2:1:atemp_diff_vel
        if(vel_row == 2)
            data_around_egg_vel_ave_filt(vel_row-1,vel_col) = nanmedian(data_around_egg_vel_ave_filt(:,vel_col));
        end
        if(abs(data_around_egg_vel_ave_filt(vel_row,vel_col)-data_around_egg_vel_ave_filt(vel_row-1,vel_col)) > 3)
            data_around_egg_vel_ave_filt(vel_row,vel_col) = data_around_egg_vel_ave_filt(vel_row-1,vel_col);
        end
    end
end



data_around_egg_prob_filt = data_around_egg_prob;

[atemp_diff_prob btemp_diff_prob] = size(data_around_egg_prob);
for vel_col = 1:1:btemp_diff_prob
    for vel_row = 2:1:atemp_diff_prob
        if(vel_row == 2)
            data_around_egg_prob_filt(vel_row-1,vel_col) = nanmedian(data_around_egg_prob_filt(:,vel_col));
        end
        if(abs(data_around_egg_prob_filt(vel_row,vel_col)-data_around_egg_prob_filt(vel_row-1,vel_col)) > 3)
            data_around_egg_prob_filt(vel_row,vel_col) = data_around_egg_prob_filt(vel_row-1,vel_col);
        end
    end
end


%% create seperate vectors for eggs on different substrates
[a b] = find(data_around_egg_sub(2400,:) == 500);
data_around_egg_500 = data_around_egg(:,b);
data_around_egg_500_sub = data_around_egg_sub(:,b);

[a b] = find(data_around_egg_sub(2400,:) == 200);
data_around_egg_200 = data_around_egg(:,b);
data_around_egg_200_sub = data_around_egg_sub(:,b);


[a b] = find(data_around_egg_sub(2400,:) == 0);
data_around_egg_0 = data_around_egg(:,b);
data_around_egg_0_sub = data_around_egg_sub(:,b);









figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar((nanmean(data_around_egg_vel_ave_filt')'), (nansem(data_around_egg_vel_ave_filt')'), time_base_egg', [.5 .5 .5])
    end
    %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
    plot(time_base_egg,(nanmean(data_around_egg_vel_ave_filt,2)),'-k');
end
ylabel('vel smoothed 1s subsampled then filtered, mm/sec');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


%data_around_egg_vel_noave_filt_sm = smooth(data_around_egg_vel_noave_filt,10);
figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar(nanmean(data_around_egg_vel_noave_filt')', nansem(data_around_egg_vel_noave_filt')', time_base_egg', [.5 .5 .5])
    end
    %plot(time_base_egg,data_around_egg_vel,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_vel_noave_filt,2),'-k');
end
ylabel('vel filtered, mm/sec');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar(nanmean(data_around_egg_prob_filt')', nansem(data_around_egg_prob_filt')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg_body,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_prob_filt,2),'-k');
end
ylabel('proboscis length filtered');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar(nanmean(data_around_egg_body_x')', nansem(data_around_egg_body_x')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg_body_x,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_body_x,2),'-k');
end
ylabel('neck to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');






figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar(nanmean(data_around_egg_angle')', nansem(data_around_egg_angle')', time_base_egg', [.5 .5 .5])
    end
    %  plot(time_base_egg,data_around_egg_angle,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_angle,2),'-k');
end
ylabel('neck to tip angle');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_body))
    if(ppp > 1)
        
        patch_errorbar(nanmean(data_around_egg')', nansem(data_around_egg')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg,2),'-k');
end
ylabel('signal around egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
end

[pp_plain,ppp_plain] = size(data_around_egg_0);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_0))
    if(ppp_plain > 1)
        
        patch_errorbar(nanmean(data_around_egg_0')', nansem(data_around_egg_0')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_0,2),'-k');
end
ylabel('signal around 0mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
end


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_0))
    if(ppp_plain > 1)
        
        patch_errorbar(nanmean(data_around_egg_0_sub')', nansem(data_around_egg_0_sub')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_0_sub,2),'-k');
end
ylabel('substrate around 0mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[0,600])
end

















[pp_plain,ppp_200] = size(data_around_egg_200);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_200))
    if(ppp_200 > 1)
        
        patch_errorbar(nanmean(data_around_egg_200')', nansem(data_around_egg_200')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_200,2),'-k');
end
ylabel('signal around 200mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
end


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_200))
    if(ppp_200 > 1)
        
        patch_errorbar(nanmean(data_around_egg_200_sub')', nansem(data_around_egg_200_sub')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_200_sub,2),'-k');
end
ylabel('substrate around 200mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[0,600])
end






[pp_plain,ppp_500] = size(data_around_egg_500);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_500))
    if(ppp_500 > 1)
        
        patch_errorbar(nanmean(data_around_egg_500')', nansem(data_around_egg_500')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_500,2),'-k');
end
ylabel('signal around 500mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
end


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_egg_500))
    if(ppp_500 > 1)
        
        patch_errorbar(nanmean(data_around_egg_500_sub')', nansem(data_around_egg_500_sub')', time_base_egg', [.5 .5 .5])
    end
    % plot(time_base_egg,data_around_egg,'color',[.75,.75,.75]);
    plot(time_base_egg,nanmean(data_around_egg_500_sub,2),'-k');
end
ylabel('substrate around 5000mM laid egg');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[0,600])
end







