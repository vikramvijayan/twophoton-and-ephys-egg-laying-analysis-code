figure; axA = subplot(2,1,1); hold on;

set(gca,'TickDir','out');

yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end),'-k')
ylabel('Vm');

yyaxis right;
ylabel('wheel position and substrate');

diff_filter_wheel = diff(recording.movie1.filtered_wheel);
colorz = 'rg';
for j = 1:5:(length(recording.movie1.filtered_wheel)-5)
    diff_wheel = (recording.movie1.filtered_wheel(j+5)-recording.movie1.filtered_wheel(j));
    if(abs(diff_wheel) < pi)
        line([recording.movie1.time_stamps(j), recording.movie1.time_stamps(j+5)], [recording.movie1.filtered_wheel(j), recording.movie1.filtered_wheel(j+5)],'LineWidth',.01,'Color',colorz(recording.movie1.sucrose(j)+1));
    end
end

scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),'k','filled');



axB = subplot(2,1,2); hold on; yyaxis right;
xlabel('Time (sec)');

plot
ylabel('body pos (AU)');

yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end),'-k');
ylabel('Vm');

linkaxes([axA, axB],'x');
%%%%%%%%%%%%%%

figure; 

ax1 = subplot(3,1,1);(recording.movie1.time_stamps, (recording.movie1.body1+recording.movie1.body2)./recording.movie1.body1,'-b')
hold on; yyaxis right;
plot(recording.movie1.time_stamps, (recording.movie1.body1+recording.movie1.body2)./recording.movie1.body1);
ylabel('body pos (AU)');
yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end))
ylabel('Vm');

ax2 =  subplot(3,1,2); hold on; yyaxis right;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH2_patch(1:1:end))
ylabel('Ch2');
yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end))
ylabel('Vm');

ax3 = subplot(3,1,3); hold on; yyaxis right;
plot(recording.abf.Time_s,recording.abf.Puff)
ylabel('puff');
yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end))
ylabel('Vm');
xlabel('Time (sec)');

linkaxes([ax1, ax2, ax3],'x');

%%%%%%%%%%%%%%

figure; 

ax1 = subplot(2,1,1);
hold on; yyaxis right;
plot(recording.movie1.time_stamps, (recording.movie1.body1+recording.movie1.body2)./recording.movie1.body1);
ylabel('body pos (AU)');
yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch(1:1:end)); ylim([-1 1])
ylabel('Vm');

ax2 =  subplot(2,1,2); hold on; yyaxis right;
plot(recording.movie1.time_stamps, (recording.movie1.body1+recording.movie1.body2)./recording.movie1.body1);
ylabel('body pos (AU)');
yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch_spikes_conv(1:1:end)); ylim([0 70])
ylabel('Vm');

linkaxes([ax1, ax2],'x');


%%%%%%%%%%%%%%

figure; axA = subplot(2,1,1); hold on;

set(gca,'TickDir','out');

yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch_spikes_conv(1:1:end),'-k'); ylim([0 70])
ylabel('Vm');

yyaxis right;
ylabel('wheel position and substrate');

diff_filter_wheel = diff(recording.movie1.filtered_wheel);
colorz = 'rg';
for j = 1:5:(length(recording.movie1.filtered_wheel)-5)
    diff_wheel = (recording.movie1.filtered_wheel(j+5)-recording.movie1.filtered_wheel(j));
    if(abs(diff_wheel) < pi)
        line([recording.movie1.time_stamps(j), recording.movie1.time_stamps(j+5)], [recording.movie1.filtered_wheel(j), recording.movie1.filtered_wheel(j+5)],'LineWidth',.01,'Color',colorz(recording.movie1.sucrose(j)+1));
    end
end

scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),'k','filled');



axB = subplot(2,1,2); hold on; yyaxis right;
xlabel('Time (sec)');

plot(recording.movie1.time_stamps, (recording.movie1.body1+recording.movie1.body2)./recording.movie1.body1,'-b')
ylabel('body pos (AU)');

yyaxis left;
plot(recording.abf.Time_s(1:1:end),recording.abf.CH1_patch_spikes_conv(1:1:end),'-k'); ylim([0 70])
ylabel('Vm');

linkaxes([axA, axB],'x');