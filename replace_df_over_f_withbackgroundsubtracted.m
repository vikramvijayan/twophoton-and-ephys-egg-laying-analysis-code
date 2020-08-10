% the way to check if this function was run is by looking at the last ROI
% (and seeing that it is all 0)
function recording = replace_df_over_f_withbackgroundsubtracted(recording)

for i =1:1:length(recording.tseries)
    for  k =1:1:(recording.tseries(i).ROI_number)
        recording.tseries(i).mean_in_ROI(k,:) = recording.tseries(i).mean_in_ROI(k,:) - recording.tseries(i).mean_in_ROI(recording.tseries(i).ROI_number,:);
        low_val = prctile(recording.tseries(i).mean_in_ROI(k,:),5);
        recording.tseries(i).df_over_f(k,:) = (recording.tseries(i).mean_in_ROI(k,:)-low_val)./low_val;
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