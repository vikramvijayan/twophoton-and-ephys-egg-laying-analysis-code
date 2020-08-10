function [transition_structure, transition_structure_notperfly] = group_plot_signal_at_transitions(recordings_list, ROI_num, filter_out_chrimson_data, filter_out_egg_times, backsub)
% ROI_num = -1 use patch
%% plot signal after transitions

% this is all transitions and never plots a second transitions
% data is stopped when there is another transition in the positive or
% negative portions of this plot

    
    transition_structure_notperfly.s200_to_p    = [];
    transition_structure_notperfly.p_to_s200    = [];
    transition_structure_notperfly.s500_to_p    = [];
    transition_structure_notperfly.p_to_s500    = [];
    transition_structure_notperfly.s500_to_s200 = [];
    transition_structure_notperfly.s200_to_s500 = [];
    
    

for rec_index = 1:1:length(recordings_list)
    %for rec_index = 1:1:12
    modifiedStr = strrep([char(recordings_list(rec_index))], '.mat', '_stripped.mat');
    %   modifiedStr = [char(recordings_list(rec_index))];
    
    recording = loadsinglerecording(modifiedStr);
    
    if(ROI_num(rec_index) ==-1)
        patch = 1;
    else
        patch = 0;
    end
    if(patch ==1)
        recording.tseries = 1;
    end
    
    [recording] = processes_more_DLC_variables(recording);
    if(backsub == 1)
        [recording] = replace_df_over_f_withbackgroundsubtracted(recording);
    end
    
    %recording.abf.CH1_patch_spikes_conv = 50000.*recording.abf.CH1_patch_spikes_conv_area_rect'./5;
    
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
                [time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num(rec_index),:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            end
            
            if(patch ==1)
                [time_base_to_return, data_to_average_interp] = average_around_event(recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
            end
            
            if(patch ==0 && filter_out_chrimson_data == 1)
                
                %tseries_no_laser = interp1(recording.abf.Time_s, recording.abf.no_laser, recording.tseries(ts).Time_s,'previous')';
                %[time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num(rec_index),:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2, recording.abf.Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
                data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
            end
            
            if(patch ==0 && filter_out_egg_times == 1)
                
                % tseries_no_egg = interp1(recording.movie1.time_stamps, recording.movie1.no_egg, recording.tseries(ts).Time_s,'previous')';
                %[time_base_to_return, data_to_average_interp] = average_around_event(recording.tseries(ts).df_over_f(ROI_num(rec_index),:), recording.tseries(ts).Time_s, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                [time_base_to_return_laser, data_to_average_interp_noegg] = average_around_event(recording.movie1.no_egg2, recording.movie1.time_stamps, recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_noegg(data_to_average_interp_noegg>0) = NaN;
                data_to_average_interp_noegg(data_to_average_interp_noegg==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_noegg;
            end
            
            if(patch ==1  && filter_out_chrimson_data == 1)
                % [time_base_to_return, data_to_average_interp] = average_around_event(recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                [time_base_to_return_laser, data_to_average_interp_laser] = average_around_event(recording.abf.no_laser2(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)), recording.movie1.time_stamps(b(i)), (startt:.01:stopt)-recording.movie1.time_stamps(b(i)));
                data_to_average_interp_laser(data_to_average_interp_laser>0) = NaN;
                data_to_average_interp_laser(data_to_average_interp_laser==0) = 1;
                data_to_average_interp = data_to_average_interp.*data_to_average_interp_laser;
                
            end
            
            
            if(patch ==1 && filter_out_egg_times == 1)
                
                disp 'Code ot yet written'
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
    
    
    for r = 1:1:rr
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) == 0) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 1 .75]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %qq = out_signal(r,:);
            % qq(out_sucrose(r,:) == 200) = NaN;
            %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            p_to_s200 = [p_to_s200, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
            
            %             qq = out_signal(r,:);
            %             qq(out_sucrose(r,:) ~= 200) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 1 .75]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %             qq = out_signal(r,:);
            %             qq(out_sucrose(r,:) ~= 0) = NaN;
            %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            s200_to_p = [s200_to_p, r];
            
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 0)
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 500) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 .75 1]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 0) = NaN;
            %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            p_to_s500 = [p_to_s500, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 0 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 500) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 .75 1]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 0) = NaN;
            %plot(-240:.1:900,qq,'color',[1 .75 .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            s500_to_p = [s500_to_p, r];
            
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 500 && recording.movie1.sucrose(store_trans(r,2)-1) == 200)
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 500) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 .75 1]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 200) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 1  .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            s200_to_s500 = [s200_to_s500, r];
        end
        
        if(recording.movie1.sucrose(store_trans(r,2)) == 200 && recording.movie1.sucrose(store_trans(r,2)-1) == 500)
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 500) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 .75 1]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            %qq = out_signal(r,:);
            %qq(out_sucrose(r,:) ~= 200) = NaN;
            %plot(-240:.1:900,qq,'color',[.75 1 .75 ]);
            %plot(-240:.1:900,qq,'color',[220 220 220]./255);
            
            s500_to_s200 = [s500_to_s200, r];
            
        end
        
    end
    
    transition_structure(rec_index).s200_to_p = nanmean(out_signal(s200_to_p,:),1)';
    transition_structure(rec_index).p_to_s200 = nanmean(out_signal(p_to_s200,:),1)';
    transition_structure(rec_index).s500_to_p = nanmean(out_signal(s500_to_p,:),1)';
    transition_structure(rec_index).p_to_s500 = nanmean(out_signal(p_to_s500,:),1)';
    transition_structure(rec_index).s500_to_s200 = nanmean(out_signal(s500_to_s200,:),1)';
    transition_structure(rec_index).s200_to_s500 = nanmean(out_signal(s200_to_s500,:),1)';
    
    transition_structure_notperfly.s200_to_p    = [transition_structure_notperfly.s200_to_p ; (out_signal(s200_to_p,:))];
    transition_structure_notperfly.p_to_s200    = [transition_structure_notperfly.p_to_s200 ;(out_signal(p_to_s200,:))];
    transition_structure_notperfly.s500_to_p    = [transition_structure_notperfly.s500_to_p ;(out_signal(s500_to_p,:))];
    transition_structure_notperfly.p_to_s500    = [transition_structure_notperfly.p_to_s500 ;(out_signal(p_to_s500,:))];
    transition_structure_notperfly.s500_to_s200 = [transition_structure_notperfly.s500_to_s200 ;(out_signal(s500_to_s200,:))];
    transition_structure_notperfly.s200_to_s500 = [transition_structure_notperfly.s200_to_s500 ;(out_signal(s200_to_s500,:))];
    
end


f_200_to_0 = figure;
f_0_to_200 = figure;
f_200_to_500 = figure;
f_500_to_200 = figure;
f_500_to_0 = figure;
f_0_to_500 = figure;


figure(f_200_to_0); hold on; title(['200 mM sucrose to plain transitions']);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_0_to_200); hold on; title(['plain to 200 mM sucrose transitions' ]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_500_to_0); hold on; title(['500 mM sucrose to plain transitions' ]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_0_to_500); hold on; title(['plain to 500 mM sucrose transitions' ]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_200_to_500); hold on; title(['200 mM sucrose to 500 mM sucrose transitions' ]);
%ylim([-0.5 3])
set(gca,'TickDir','out');

figure(f_500_to_200); hold on; title(['500 mM sucrose to 200 mM sucrose transitions' ]);
%ylim([-0.5 3])
set(gca,'TickDir','out');





out_signal = [transition_structure(1:end).s200_to_p];

figure(f_200_to_0); hold on;
if(~isempty(out_signal))
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[33 113 181]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[107 174 214]./255,'linewidth',2);
    
    xlim([-45,45]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = [transition_structure(1:end).p_to_s200];

if(~isempty(out_signal))
    
    figure(f_0_to_200); hold on;
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[107 174 214]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[33 113 181]./255,'linewidth',2);
    
    xlim([-45,45]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = [transition_structure(1:end).s500_to_p];
if(~isempty(out_signal))
    
    figure(f_500_to_0); hold on;
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[8 48 107]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[107 174 214]./255,'linewidth',2);
    
    xlim([-45,45]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = [transition_structure(1:end).p_to_s500];
if(~isempty(out_signal))
    
    figure(f_0_to_500); hold on;
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[107 174 214]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[8 48 107]./255,'linewidth',2);
    
    xlim([-45,45]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = [transition_structure(1:end).s500_to_s200];
if(~isempty(out_signal))
    
    figure(f_500_to_200); hold on;
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[8 48 107]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[33 113 181]./255,'linewidth',2);
    
    xlim([-45,45]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = [transition_structure(1:end).s200_to_s500];
if(~isempty(out_signal))
    
    figure(f_200_to_500); hold on;
    qq = nanmean(out_signal,2);
    plot(-240:.1:900,(out_signal),'color',[220 220 220]./255,'linewidth',1);
    plot(-240:.1:0,qq(1:2401),'color',[33 113 181]./255,'linewidth',2);
    plot(0:.1:900,qq((end-9000):end),'color',[8 48 107]./255,'linewidth',2);
    
    xlim([-45,45]);
end

figure(f_200_to_0);  set(gca,'ylim',[-.5,1]); pause(1); snapnow;

figure(f_0_to_200); set(gca,'ylim',[-.5,1]);  pause(1);snapnow;

figure(f_500_to_0); set(gca,'ylim',[-.5,1]);  pause(1); snapnow;

figure(f_0_to_500); set(gca,'ylim',[-.5,1]); pause(1); snapnow;

figure(f_200_to_500); set(gca,'ylim',[-.5,1]); pause(1); snapnow;

figure(f_500_to_200); set(gca,'ylim',[-.5,1]); pause(1); snapnow;


figure(f_200_to_0);  set(gca,'ylim',[0,.5]); pause(1); snapnow;

figure(f_0_to_200); set(gca,'ylim',[0,.5]); pause(1);snapnow;

figure(f_500_to_0);  set(gca,'ylim',[0,.5]);  pause(1); snapnow;

figure(f_0_to_500);  set(gca,'ylim',[0,.5]); pause(1); snapnow;

figure(f_200_to_500);  set(gca,'ylim',[0,.5]); pause(1); snapnow;

figure(f_500_to_200);  set(gca,'ylim',[0,.5]); pause(1); snapnow;


figure(f_200_to_0); xlim([-120,120]);   pause(1); snapnow;
figure(f_0_to_200); xlim([-120,120]);  pause(1); snapnow;
figure(f_500_to_0); xlim([-120,120]);  pause(1);  snapnow;
figure(f_0_to_500); xlim([-120,120]);  pause(1); snapnow;
figure(f_200_to_500); xlim([-120,120]); pause(1);  snapnow;
figure(f_500_to_200); xlim([-120,120]);   pause(1); snapnow;

out=1;

end

