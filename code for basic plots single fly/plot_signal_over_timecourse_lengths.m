% plots all lengths (4 subplots total)

function [out] = plot_signal_over_timecourse_lengths(recording, ROI_num)


fig = figure; 
set(fig,'defaultAxesColorOrder',[0,0,1 ; 0,0,0]);
axA = subplot(4,1,1);
set(gca,'FontSize',7)

hold on;
box on;
set(gca,'TickDir','out');
set(gca,'ylim',[0, 2*pi]);
title([recording.abfname ',ROI ' num2str(ROI_num)],'Interpreter','none');
%set(gca,'xticklabel',{[]});
ylabel('wheel position');
set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])

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

scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),'m','filled');

for i =1:1:length(recording.movie1.eggs)
    h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color','m' );
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
yyaxis right;
ylabel('df over f');

for i = 1:1:length(recording.tseries)
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')
end


yyaxis left;
scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),'m','filled');

axB = subplot(4,1,2);
set(gca,'FontSize',7)

hold on;
box on;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');
  

 yyaxis right;
ylabel('df over f');


for i = 1:1:length(recording.tseries)
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color',[.3,.3,.3])
end


 yyaxis left;
ylabel('len - L3 to tip PW');

plot(recording.movie1.time_stamps,recording.movie1.abd_path_length,'-','color','b')

set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])



axC = subplot(4,1,3);
set(gca,'FontSize',7)

hold on;
box on;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');
  


 yyaxis right;
ylabel('df over f');

for i = 1:1:length(recording.tseries)
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color',[.3,.3,.3])
end

 yyaxis left;
ylabel('len - L3 to tip');

plot(recording.movie1.time_stamps,recording.movie1.abd_only_length,'-','color','b')

set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])




axD = subplot(4,1,4);
set(gca,'FontSize',7)

hold on;
box on;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');
  


 yyaxis right;
ylabel('df over f');

for i = 1:1:length(recording.tseries)
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color',[.3,.3,.3])
end

 yyaxis left;
ylabel('len - N to tip');

plot(recording.movie1.time_stamps,recording.movie1.abd_length,'-','color','b')

set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])




set(gcf, 'Position', [0, 100, 1800, 150]);

%set(gcf, 'Position', [100, 200, 1400, 500]);
Gx = gcf;
Gx.Position(3:4) = Gx.Position(3:4)*4;
% Ax = gca;
% Ax.FontSize = Ax.FontSize *3;

%snapnow
% for lenplot = min(recording.abf.Time_s):1800:(max(recording.abf.Time_s))
%     set(gca,'xlim',[lenplot, lenplot+1800])
%     snapnow
% end
linkaxes([axA,axB,axC,axD],'x');

snapnow
%close all;





out=1;

end

