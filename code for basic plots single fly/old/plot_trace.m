figure; hold on;

diff_filter_wheel = diff(recording.movie1.filtered_wheel);

colorz = 'rg';
for j = 1:50:(length(recording.movie1.filtered_wheel)-50)
    diff_wheel = (recording.movie1.filtered_wheel(j+50)-recording.movie1.filtered_wheel(j));
    if(abs(diff_wheel) < pi)
        line([recording.movie1.time_stamps(j), recording.movie1.time_stamps(j+50)], [recording.movie1.filtered_wheel(j), recording.movie1.filtered_wheel(j+50)],'LineWidth',.01,'Color',colorz(recording.movie1.sucrose(j)+1));
    end
end

scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.eggs),'k','filled');
scatter(recording.movie1.time_stamps(recording.movie1.eggs), recording.movie1.filtered_wheel(recording.movie1.heavy_filt_transition),'b','filled');




for i = 1:1:length(recording.tseries)
    
    yyaxis right;
    
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f_gauss(1,:)),[.3,.3,.9])
    plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f_gauss(2,:)),[.5,.5,.9])
    
    %   plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(1,:)),'-b')
    % plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f(2,:)),'-m')
    
    %         plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f_gauss(1,:)),'-b')
    %     plot(recording.tseries(i).Time_s,(recording.tseries(i).df_over_f_gauss(2,:)),'-m')
    
    set(gca,'TickDir','out');
    
end