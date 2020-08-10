% if ROI_num = -1, use patch
% this is the most useful plotting function so far
% very streamlined

% 1- signal and wheel
% 2- abd length x
% 3- abd angle
function [out] = plot_signal_over_timecourse_streamlined_median_chrimson_testpat(recording, ROI_num)



% process more of the DLC movie
% this was predone for the patching files
% [recording] = processes_more_DLC_variables(recording);


fig = figure;
set(fig,'defaultAxesColorOrder',[0,0,1 ; 0,0,0]);
axA = subplot(4,1,1);
set(gca,'FontSize',7)

hold on;
box off;
set(gca,'TickDir','out');
%set(gca,'ylim',[0, 2*pi]);
title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
set(gca,'xticklabel',{[]});
ylabel('wheel position');
set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])
set(gca,'ylim',[0, 2*pi])

diff_filter_wheel = diff(recording.movie1.filtered_wheel);

colorz = [ 107./255,174./255,214./255; [33 113 181]./255; [33 113 181]./255; 8./255,48./255,107./255; 8./255,48./255,107./255; 8./255,48./255,107./255];
for j = 1:10:(length(recording.movie1.filtered_wheel)-10)
    diff_wheel = (recording.movie1.filtered_wheel(j+10)-recording.movie1.filtered_wheel(j));
    if(abs(diff_wheel) < pi)
        h = line([recording.movie1.time_stamps(j), recording.movie1.time_stamps(j+10)], [recording.movie1.filtered_wheel(j), recording.movie1.filtered_wheel(j+10)],'LineWidth',1,'Color',colorz(recording.movie1.sucrose(j)./100+1,:));
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end


ds = diff(recording.movie1.sucrose);
ds = [0,ds];
[a b] = find(ds);
b = [1,b, length(recording.movie1.sucrose)];

for i =2:1:(length(b)-1)
    % time between transition (i-1) and i+1 has to be greater than 4 seconds and has to be in transition zone (removes errant)
    if( (b(i+1)-b(i-1)) > 100)
        h = line([recording.movie1.time_stamps(b(i)), recording.movie1.time_stamps(b(i))], [0, 2*pi],'color',[150 150 150]./255);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

% put the dot onthe wheel position
scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),4,'m','filled');


for i =1:1:length(recording.movie1.eggs)
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[150 150 150]./255 );
    h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color','m' );
    
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end


%% more behavior if it is there
data_around_egg_behavior = [10381,10956,11456,11556,12156,19978,20044,20127;36644,36912,37535,37765,37995,41373,41419,41600;21086,21537,22112,22514,22946,26350,26398,26436;27748,27854,28189,28573,28573,29459,29908,30106];
for i =1:1:4
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[150 150 150]./255 );
    h = line([recording.movie1.time_stamps(data_around_egg_behavior(i,6)), recording.movie1.time_stamps(data_around_egg_behavior(i,6))], [0, 2*pi],'color',[128,0,128]./255);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,6)), 0,6,[128,0,128]./255,'filled');
    
end

for i =1:1:4
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[216,191,216]./255 );
    h = line([recording.movie1.time_stamps(data_around_egg_behavior(i,1)), recording.movie1.time_stamps(data_around_egg_behavior(i,1))], [0, 2*pi],'color',[216,191,216]./255);
    scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,1)), 0,6,[216,191,216]./255,'filled');
    
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end




%plot(recording.abf.Time_s, recording.abf.PWMlaser,'r');
%plot(recording.abf.Time_s, recording.abf.PWMlaser_uWpermm2,'r');

yyaxis right;

if (ROI_num > 0)
    ylabel('df over f');
    
    for i = 1:1:length(recording.tseries)
        plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')
    end
else
    ylabel('5 s rect spikes per sec');
    %plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



window = 50000; %5 sec

y = filter((1/window)*ones(1,window),1,recording.abf.CH1_patch_spikes');
y = [ y((window/2)+1:end), NaN(1,window/2)];
y = y';
    plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*y(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



end

% yyaxis left;
% % put the dot onthe wheel position
% scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),4,'m','filled');


axB = subplot(4,1,2); hold on;

set(gca,'FontSize',7)

box off;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');

set(gca,'xticklabel',{[]});
yyaxis right; 

yyaxis right; hold on;

if (ROI_num > 0)
    ylabel('df over f');
    
    for i = 1:1:length(recording.tseries)
        plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')
    end
else
    ylabel('5 s gaussian spikes per sec');
    plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')


end

axC = subplot(4,1,3); hold on;

set(gca,'FontSize',7)

box off;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');

set(gca,'xticklabel',{[]});
yyaxis right;  hold on;


if (ROI_num > 0)
    ylabel('df over f');
    
    for i = 1:1:length(recording.tseries)
        plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')
    end
else
    ylabel('10 s gaussian spikes per sec');
    %plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



window = 100000; %2 sec

w = gausswin(window); 
y = filter(w./(sum(w)),1,recording.abf.CH1_patch_spikes');

y = [ y((window/2)+1:end), NaN(1,window/2)];
y=y';

plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*y(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



end



axD = subplot(4,1,4); hold on;

set(gca,'FontSize',7)

box off;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');

set(gca,'xticklabel',{[]});
yyaxis right;  hold on;


if (ROI_num > 0)
    ylabel('df over f');
    
    for i = 1:1:length(recording.tseries)
        plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')
    end
else
    ylabel('2 s gausswin spikes per sec');
    %plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



window = 20000; %2 sec

w = gausswin(window); 
y = filter(w./(sum(w)),1,recording.abf.CH1_patch_spikes');

y = [ y((window/2)+1:end), NaN(1,window/2)];
y = y';
    plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*y(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')



end

Gx = gcf;
Gx.Position(3:4) = Gx.Position(3:4)*6;
% Ax = gca;
% Ax.FontSize = Ax.FontSize *3;
linkaxes([axA,axB,axC,axD],'x');
out=1;

end

