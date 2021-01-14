% the way to check if this function was run is by looking at the last ROI
% (and seeing that it is all 0)
function recording = replace_df_over_f_withbackgroundsubtracted_runningmean(recording)

for i =1:1:length(recording.tseries)
    for  k =1:1:(recording.tseries(i).ROI_number)
        recording.tseries(i).mean_in_ROI(k,:) = recording.tseries(i).mean_in_ROI(k,:) - recording.tseries(i).mean_in_ROI(recording.tseries(i).ROI_number,:);
        
        
        interpsignal = interp1(recording.tseries(i).Time_s,recording.tseries(i).mean_in_ROI(k,:), recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),'previous');
        interpsignal_sm = smoothdata(interpsignal,'movmean',10*60*20,'omitnan');
        
        % uses a 20 min average normalization
        tmp_var = (interpsignal  -interpsignal_sm) ./ interpsignal_sm;
        interpsignal_back = interp1(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),tmp_var, recording.tseries(i).Time_s,'previous');
        recording.tseries(i).df_over_f(k,:) = interpsignal_back;
        
        % uses a full recording normalization
          tmp_var = (interpsignal  - nanmean(interpsignal)) ./ nanmean(interpsignal);
        interpsignal_back = interp1(recording.tseries(i).Time_s(1):.1:recording.tseries(i).Time_s(end),tmp_var, recording.tseries(i).Time_s,'previous');
        recording.tseries(i).df_over_f_nosm(k,:) = interpsignal_back;

    end
    tmp = recording.tseries(i).df_over_f(recording.tseries(i).ROI_number,:);
    recording.tseries(i).df_over_f(recording.tseries(i).ROI_number,isnan(tmp)) = 0;
end



% for i =1:1:length(recording.tseries)
%     
%     recording.tseries(i).mean_in_ROI(1,:) = (recording.tseries(i).mean_in_ROI(1,:) + recording.tseries(i).mean_in_ROI(2,:))./2 - recording.tseries(i).mean_in_ROI(recording.tseries(i).ROI_number,:);
%     low_val = prctile(recording.tseries(i).mean_in_ROI(1,:),5);
%     recording.tseries(i).df_over_f(1,:) = (recording.tseries(i).mean_in_ROI(1,:)-low_val)./low_val;
%     
%     tmp = recording.tseries(i).df_over_f(recording.tseries(i).ROI_number,:);
%     recording.tseries(i).df_over_f(recording.tseries(i).ROI_number,isnan(tmp)) = 0;
% end