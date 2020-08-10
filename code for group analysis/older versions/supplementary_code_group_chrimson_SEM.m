pulse_structure = pulse_structure_median;


patch = 0;

% plotting pulse triggered averages
data_around_pulse_vel = [pulse_structure(1:end).data_around_pulse_vel];
data_around_pulse_vel_noave = [pulse_structure(1:end).data_around_pulse_vel_noave];


data_around_pulse_body = [pulse_structure(1:end).data_around_pulse_body];
data_around_pulse_body_x = [pulse_structure(1:end).data_around_pulse_body_x];
data_around_pulse_body_y = [pulse_structure(1:end).data_around_pulse_body_y];

data_around_pulse_path = [pulse_structure(1:end).data_around_pulse_path];
data_around_pulse_body_only = [pulse_structure(1:end).data_around_pulse_body_only];
data_around_pulse_body_only_x = [pulse_structure(1:end).data_around_pulse_body_only_x];
data_around_pulse_body_only_y = [pulse_structure(1:end).data_around_pulse_body_only_y];

data_around_pulse_angle = [pulse_structure(1:end).data_around_pulse_angle];
data_around_pulse_angle_only = [pulse_structure(1:end).data_around_pulse_angle_only];
data_around_pulse = [pulse_structure(1:end).data_around_pulse];

data_around_pulse_prob = [pulse_structure(1:end).data_around_pulse_prob];
data_around_pulse_prob_x = [pulse_structure(1:end).data_around_pulse_prob_x];
data_around_pulse_prob_y = [pulse_structure(1:end).data_around_pulse_prob_y];

data_around_pulse_temp = [pulse_structure(1:end).data_around_pulse_temp];
data_around_pulse_PWMlaseruW = [pulse_structure(1:end).data_around_pulse_PWMlaseruW];
data_around_pulse_puff = [pulse_structure(1:end).data_around_pulse_puff];

data_around_pulse_egglaid = [pulse_structure(1:end).data_around_pulse_egglaid];

% remove double plotting of behavior
if(~isempty(data_around_pulse))
[a b] = find(~isnan(data_around_pulse(2400,:)));

data_around_pulse_vel = data_around_pulse_vel(:,b);
data_around_pulse_vel_noave = data_around_pulse_vel_noave(:,b);

data_around_pulse_body = data_around_pulse_body(:,b);
data_around_pulse_body_x = data_around_pulse_body_x(:,b);
data_around_pulse_body_y = data_around_pulse_body_y(:,b);
data_around_pulse_path = data_around_pulse_path(:,b);
data_around_pulse_body_only = data_around_pulse_body_only(:,b);
data_around_pulse_body_only_x = data_around_pulse_body_only_x(:,b);
data_around_pulse_body_only_y = data_around_pulse_body_only_y(:,b);
data_around_pulse_angle = data_around_pulse_angle(:,b);
data_around_pulse_angle_only = data_around_pulse_angle_only(:,b);
data_around_pulse = data_around_pulse(:,b);
data_around_pulse_prob = data_around_pulse_prob(:,b);
data_around_pulse_prob_x = data_around_pulse_prob_x(:,b);
data_around_pulse_prob_y = data_around_pulse_prob_y(:,b);
data_around_pulse_temp = data_around_pulse_temp(:,b);
data_around_pulse_PWMlaseruW = data_around_pulse_PWMlaseruW(:,b);
data_around_pulse_puff = data_around_pulse_puff(:,b);

data_around_pulse_egglaid = data_around_pulse_egglaid(:,b);
end

[pp,ppp] = size(data_around_pulse_body);

time_base_pulse = [-240:.1:240];

% flter out large jumps in velocity and proboscis
data_around_pulse_vel_noave_filt = data_around_pulse_vel_noave;

[atemp_diff_vel btemp_diff_vel] = size(data_around_pulse_vel_noave);
for vel_col = 1:1:btemp_diff_vel
    for vel_row = 2:1:atemp_diff_vel
        if(vel_row == 2)
            data_around_pulse_vel_noave_filt(vel_row-1,vel_col) = nanmedian(data_around_pulse_vel_noave_filt(:,vel_col));
        end
        if(abs(data_around_pulse_vel_noave_filt(vel_row,vel_col)-data_around_pulse_vel_noave_filt(vel_row-1,vel_col)) > 1)
            data_around_pulse_vel_noave_filt(vel_row,vel_col) = data_around_pulse_vel_noave_filt(vel_row-1,vel_col);
        end
    end
end

data_around_pulse_vel_ave_filt = data_around_pulse_vel;

[atemp_diff_vel btemp_diff_vel] = size(data_around_pulse_vel);
for vel_col = 1:1:btemp_diff_vel
    for vel_row = 2:1:atemp_diff_vel
        if(vel_row == 2)
            data_around_pulse_vel_ave_filt(vel_row-1,vel_col) = nanmedian(data_around_pulse_vel_ave_filt(:,vel_col));
        end
        if(abs(data_around_pulse_vel_ave_filt(vel_row,vel_col)-data_around_pulse_vel_ave_filt(vel_row-1,vel_col)) > 1)
            data_around_pulse_vel_ave_filt(vel_row,vel_col) = data_around_pulse_vel_ave_filt(vel_row-1,vel_col);
        end
    end
end

data_around_pulse_prob_filt = data_around_pulse_prob;

[atemp_diff_prob btemp_diff_prob] = size(data_around_pulse_prob);
for vel_col = 1:1:btemp_diff_prob
    for vel_row = 2:1:atemp_diff_prob
        if(vel_row == 2)
            data_around_pulse_prob_filt(vel_row-1,vel_col) = nanmedian(data_around_pulse_prob_filt(:,vel_col));
        end
        if(abs(data_around_pulse_prob_filt(vel_row,vel_col)-data_around_pulse_prob_filt(vel_row-1,vel_col)) > 3)
            data_around_pulse_prob_filt(vel_row,vel_col) = data_around_pulse_prob_filt(vel_row-1,vel_col);
        end
    end
end



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
    if(ppp > 1)
    patch_errorbar(nanmean(data_around_pulse_vel_noave_filt')', nansem(data_around_pulse_vel_noave_filt')', time_base_pulse', [.5 .5 .5]);
    end
    %plot(time_base_pulse,data_around_pulse_vel,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_vel_noave_filt,2),'-k');
end
ylabel('vel filtered, mm/sec');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
    if(ppp > 1)
        
        patch_errorbar((nanmean(data_around_pulse_vel_ave_filt')'), (nansem(data_around_pulse_vel_ave_filt')'), time_base_pulse', [.5 .5 .5])
    end
    %plot(time_base_pulse,data_around_pulse_vel,'color',[.75,.75,.75]);
    plot(time_base_pulse,(nanmean(data_around_pulse_vel_ave_filt,2)),'-k');
end
ylabel('vel smoothed 1s subsampled then filtered, mm/sec');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');





figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_body_x')', nansem(data_around_pulse_body_x')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse_body_x,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_body_x,2),'-k');
end
ylabel('neck to tip length, X only');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);



figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_angle')', nansem(data_around_pulse_angle')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_angle,2),'-k');
end
ylabel('neck to tip angle');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);

figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse_prob_filt')', nansem(data_around_pulse_prob_filt')', time_base_pulse', [.5 .5 .5])
        end
    %  plot(time_base_pulse,data_around_pulse_angle,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse_prob_filt,2),'-k');
end
ylabel('proboscis length filtered');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
snapnow;
set(gca,'xlim',[-60,60]);


figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'TickDir','out');
if(~isempty(data_around_pulse_body))
        if(ppp > 1)

    patch_errorbar(nanmean(data_around_pulse')', nansem(data_around_pulse')', time_base_pulse', [.5 .5 .5])
        end
    % plot(time_base_pulse,data_around_pulse,'color',[.75,.75,.75]);
    plot(time_base_pulse,nanmean(data_around_pulse,2),'-k');
end
ylabel('signal around pulse');
xlabel('time (sec)'); %since peak is neg, patch comes first. however...
axis manual;
line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
if(~patch)
    set(gca,'ylim',[-.5,1.5])
    %set(gca,'ylim',[-.5,3.5])
end
snapnow;
set(gca,'xlim',[-60,60]);





% 
% figure; hold on; %title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
 figure; plot(time_base_pulse,10.*smooth(1.*nansum(data_around_pulse_egglaid,2),50))




% 
% set(gca,'TickDir','out');
% if(~isempty(data_around_pulse_body))
%         if(ppp > 1)
% 
%     patch_errorbar(nanmean(data_around_pulse_egglaid')', nansem(data_around_pulse_egglaid')', time_base_pulse', [.5 .5 .5])
%         end
%     % plot(time_base_pulse,data_around_pulse,'color',[.75,.75,.75]);
%     plot(time_base_pulse,nanmean(data_around_pulse_egglaid,2),'-k');
% end
% ylabel('signal around pulse egglaid');
% xlabel('time (sec)'); %since peak is neg, patch comes first. however...
% axis manual;
% line('Xdata',[0,0],'YData',[-1000,1000],'Color','c');
% if(~patch)
%     set(gca,'ylim',[0,1])
%     %set(gca,'ylim',[-.5,3.5])
% end
% snapnow;
%set(gca,'xlim',[-60,60]);