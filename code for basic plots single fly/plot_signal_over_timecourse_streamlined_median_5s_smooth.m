% if ROI_num = -1, use patch
% this is the most useful plotting function so far
% very streamlined

% 1- signal and wheel
% 2- abd length x
% 3- abd angle
function [out] = plot_signal_over_timecourse_streamlined_median_5s_smooth(recording, ROI_num)



% process more of the DLC movie
% this was predone for the patching files
[recording] = processes_more_DLC_variables(recording);


fig = figure; 
set(fig,'defaultAxesColorOrder',[0,0,1 ; 0,0,0]);
axA = subplot(5,1,1);
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

% for i =2:1:(length(b)-1)
%     % time between transition (i-1) and i+1 has to be greater than 4 seconds and has to be in transition zone (removes errant)
%     if( (b(i+1)-b(i-1)) > 100)
%         h = line([recording.movie1.time_stamps(b(i)), recording.movie1.time_stamps(b(i))], [0, 2*pi],'color',[150 150 150]./255);
%         set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     end
% end

% put the dot onthe wheel position
%scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),4,'m','filled');

axD = subplot(5,1,2);
hold on;
set(gca,'FontSize',7)
yyaxis right;

for i =1:1:length(recording.movie1.eggs)
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[150 150 150]./255 );
    h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color','r' );
    scatter(recording.movie1.time_stamps(recording.movie1.eggs(i)), 0,6,'r','filled');

    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
% 

%% more behavior if it is there
 %data_around_egg_behavior = [48532,49640,49910,50617,50623,52586,52830,52830;59818,61572,61700,62519,62801,65160,65370,65370];


 % 2019_08_01_0003
 %data_around_egg_behavior = [16075,17948,18348,20082,20270,20757,20761,20944,0,0,0,0;44292,44980,45800,46730,47075,47901,47982,48156,0,0,0,0;54738,55232,55736,56670,56760,57714,57773,57900,0,0,0,0;72718,73151,73655,74164,74168,75034,75560,75670,0,0,0,0;64070,64685,65233,66042,66054,66470,66963,67193,0,0,0,0;96344,96857,97318,98103,98911,99425,99769,99945,0,0,0,0;102472,103100,103613,104478,104605,105598,105612,105879,0,0,0,0;108118,108464,109095,109886,110134,110619,110619,110893,0,0,0,0;112893,113110,113638,114327,114693,116148,116148,116299,0,0,0,0;118694,118884,119896,120790,120877,122120,122120,122276,0,0,0,0;139915,140472,141146,141644,142421,143266,143266,143390,0,0,0,0];

 %2019_06_19_0003
 %data_around_egg_behavior = [29242,29560,30086,30770,30770,31805,32093,32341,0,0,0,0;44415,44956,45383,45800,45800,47360,47455,47606,0,0,0,0;76005,79248,79917,80413,80431,82633,82682,82727,0,0,0,0;159234,159248,160253,160979,161591,162511,162622,162709,0,0,0,0];
 
 % chrimson eggs 202007170007
%data_around_egg_behavior = [29689,29689,29689,29689,32152,33409,33820,29689; 55209,55209,55209,55209,57044,58729,59214,55209];


%2019 08050001
%data_around_egg_behavior = [60679,61861,62145,62631,63489,64113,65028,65183,0,0,0,0;115236,116185,116346,116801,117401,118145,118841,118984,0,0,0,0;131845,132556,132854,133679,133857,134207,134850,134950,0,0,0,0;153443,154016,154475,154933,155292,156209,156209,156407,0,0,0,0;138655,140924,141268,141842,141842,142529,143170,143344,0,0,0,0;163073,163348,163508,164219,164907,165526,165596,165820,0,0,0,0;182401,183341,183501,184409,184409,185403,185444,185600,0,0,0,0];

%2019 08060003
%data_around_egg_behavior =  [51864,52410,52780,53478,53478,54656,55156,55274,0,0,0,0;163440,163600,164075,164740,164740,166112,168780,168855,0,0,0,0;178923,179566,179947,180268,180268,181939,181939,182128,0,0,0,0];


%2019 08060004
%data_around_egg_behavior = [13851,14548,14896,16546,16546,17083,17083,17240,0,0,0,0;29412,29732,30125,30826,30826,32448,32448,32612,0,0,0,0;41326,41629,41995,42563,42563,44244,44244,44380,0,0,0,0;54115,54334,54906,55394,55394,56974,56974,57225,0,0,0,0;65676,66130,66584,67060,67060,70483,70483,70656,0,0,0,0];
   
%2020 10150003 patch

%data_around_egg_behavior = [33840,34356,34919,36295,36649,38872,38872,39025];

%2020 10150001 patch
%data_around_egg_behavior = [90029,90837,91316,92326,92831,97245,97245,97339];


data_around_egg_behavior = [];

[num_eggs ~] = size(data_around_egg_behavior);

for i =1:1:num_eggs
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[150 150 150]./255 );
    h = line([recording.movie1.time_stamps(data_around_egg_behavior(i,6)), recording.movie1.time_stamps(data_around_egg_behavior(i,6))], [0, 2*pi],'color',[128,0,128]./255);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,6)), 0,6,[128,0,128]./255,'filled');
    
end

for i =1:1:num_eggs
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[216,191,216]./255 );
    h = line([recording.movie1.time_stamps(data_around_egg_behavior(i,1)), recording.movie1.time_stamps(data_around_egg_behavior(i,1))], [0, 2*pi],'color',[216,191,216]./255);
    scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,1)), 0,6,[216,191,216]./255,'filled');
    
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

for i =1:1:num_eggs
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[216,191,216]./255 );
    h = line([recording.movie1.time_stamps(data_around_egg_behavior(i,5)), recording.movie1.time_stamps(data_around_egg_behavior(i,5))], [0, 2*pi],'color',[216,191,216]./255);
    scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,5)), 0,6,[240,107,168]./255,'filled');
    
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end


%more_bends = [26319,27644,35152,36983,40419];
more_bends = [];

for i =1:1:length(more_bends)
    %h = line([recording.movie1.time_stamps(recording.movie1.eggs(i)), recording.movie1.time_stamps(recording.movie1.eggs(i))], [0, 2*pi],'color',[216,191,216]./255 );
    h = line([recording.movie1.time_stamps(more_bends(i)), recording.movie1.time_stamps(more_bends(i))], [0, 2*pi],'color',[128,0,128]./255);
    %scatter(recording.movie1.time_stamps(data_around_egg_behavior(i,1)), 0,6,[216,191,216]./255,'filled');
    
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

%plot(recording.abf.Time_s, recording.abf.PWMlaser,'r');
plot(recording.abf.Time_s, recording.abf.PWMlaser_uWpermm2./100,'color',[255,165,0]./255);

yyaxis left;

if (ROI_num > 0)
    ylabel('df over f mean normalized 20 min window back sub, 2s smooth');
    %ylabel('df over f fo normalized, 5s smooth');
    %ylabel('mean pixel value in ROI no normalizing, 2s smooth');
    %ylabel('df over f, fo is mean, 1s smooth');

    for i = 1:1:length(recording.tseries)             
        
        %interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).df_over_f(ROI_num,:)./nanmean(recording.tseries(i).df_over_f(ROI_num,:)), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        % interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).df_over_f(ROI_num,:), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        %interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).mean_in_ROI(ROI_num,:), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        %interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).mean_in_ROI(ROI_num,:)./nanmean(recording.tseries(i).mean_in_ROI(ROI_num,:)), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        
        % interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).mean_in_ROI(ROI_num,:)./nanmean(recording.tseries(i).mean_in_ROI(ROI_num,:))-1, recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        
        %plot(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),smoothdata(interpsignal,'movmean',50,'omitnan'),'-','color','b')
        
        
        % interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).df_over_f(ROI_num,:), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        % interpsignal_sm =          interpsignal./smoothdata(interpsignal,'movmean',12000,'omitnan');
        % plot(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),smoothdata(interpsignal_sm,'movmean',20,'omitnan'),'-','color','k')
        
       % plot(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),smoothdata(interpsignal,'movmean',20,'omitnan'),'-','color','k')

       
               interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).df_over_f(ROI_num,:), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
                plot(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),smoothdata(interpsignal,'movmean',20,'omitnan'),'-','color','k')

    end
else
    %ylabel('spikes per sec');
    ylabel('corrected Vm');

    %plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')
    plot(recording.abf.Time_s(recording.time_to_use(1)*10000:10:floor(recording.time_to_use(2)*10000)),recording.abf.CH1_patch(recording.time_to_use(1)*10000:10:floor(recording.time_to_use(2)*10000))-13,'-','color','k')

end

  yyaxis left;
% put the dot onthe wheel position
%scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),4,'m','filled');

axE = subplot(5,1,3);

hold on;
box off;
set(gca,'FontSize',7)
set(gca,'TickDir','out');

if(ROI_num < 0)
        ylabel('spikes per sec');
    %ylabel('corrected Vm');

    plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),10000.*recording.abf.CH1_patch_spikes_conv_area_rect(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')
    %plot(recording.abf.Time_s(recording.time_to_use(1)*10000:10:floor(recording.time_to_use(2)*10000)),recording.abf.CH1_patch(recording.time_to_use(1)*10000:10:floor(recording.time_to_use(2)*10000))-13,'-','color','k')

end


axB = subplot(5,1,4);
set(gca,'FontSize',7)

hold on;
box off;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');
  
set(gca,'xticklabel',{[]});


%  yyaxis right;
% ylabel('df over f');
% 
% if (ROI_num > 0)
%     for i = 1:1:length(recording.tseries)
%         plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')      
%     end
% else
%     plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')
% end


 yyaxis left;
 
%ylabel('abd length - N to tip X / median smoothed');
%plot(recording.movie1.time_stamps,smooth(recording.movie1.abd_x_neck_tip./(nanmedian(recording.movie1.abd_x_neck_tip)),25),'-','color','b');
ylabel('abd length - N to tip X / median 2s smooth');
% sm = smoothdata((recording.movie1.abd_x_neck_tip./(nanmedian(recording.movie1.abd_x_neck_tip))),'movmean',125,'omitnan');
% plot(recording.movie1.time_stamps,sm,'-','color','b')
sm = smoothdata((recording.movie1.abd_x_neck_tip./(nanmedian(recording.movie1.abd_x_neck_tip))),'movmean',50,'omitnan');
plot(recording.movie1.time_stamps,sm,'-','color','k')


set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])






axC = subplot(5,1,5);
set(gca,'FontSize',7)

hold on;
box off;
set(gca,'TickDir','out');
% set(gca,'ylim',[-.5, 3]);
% ylabel('ROI 1');
  


%  yyaxis right;
% ylabel('df over f');
% 
% if (ROI_num > 0)
%     for i = 1:1:length(recording.tseries)
%         plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(ROI_num,:)),'-','color','k')      
%     end
% else
%     plot(recording.abf.Time_s(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),recording.abf.CH1_patch_spikes_conv(recording.time_to_use(1)*10000:100:floor(recording.time_to_use(2)*10000)),'-','color','k')
% end


 yyaxis left;
%ylabel('abd angle smoothed');
%plot(recording.movie1.time_stamps,smooth(recording.movie1.abd_angle,25),'-','color','b')

ylabel('abd angle, 2s smooth');
% sm = smoothdata(recording.movie1.abd_angle,'movmean',125,'omitnan');
% plot(recording.movie1.time_stamps,sm,'-','color','b')
sm = smoothdata((360/(2*pi)).*recording.movie1.abd_angle,'movmean',50,'omitnan');
plot(recording.movie1.time_stamps,sm,'-','color','k')

%set(gca,'ylim',[.5, pi/2])

set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])

set(gca,'xlim',[min(recording.abf.Time_s), max(recording.abf.Time_s)])
%set(gca,'ylim',[.5, 2*pi/3])


set(gcf, 'Position', [0, 100, 1800, 150]);

%set(gcf, 'Position', [100, 200, 1400, 500]);
Gx = gcf;
Gx.Position(3:4) = Gx.Position(3:4)*6;
% Ax = gca;
% Ax.FontSize = Ax.FontSize *3;
linkaxes([axA,axB,axC, axD, axE],'x');


% for egg_cnt = 1:1:length(recording.movie1.eggs)
%     qq = max(recording.movie1.eggs(egg_cnt)-600*25,1);
%     qq2 = min(recording.movie1.eggs(egg_cnt)+240*25,length(recording.movie1.time_stamps));
%     axis(axF);
%     xlim([recording.movie1.time_stamps(qq),recording.movie1.time_stamps(qq2)  ]);
%     snapnow
% end
% for lenplot = min(recording.abf.Time_s):1800:(max(recording.abf.Time_s))
%     set(gca,'xlim',[lenplot, lenplot+1800])
%     snapnow
% end


%close all;



out=1;

end

