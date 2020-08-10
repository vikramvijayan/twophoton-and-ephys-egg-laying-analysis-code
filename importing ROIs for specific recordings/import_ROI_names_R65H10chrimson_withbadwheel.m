% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_R65H10chrimson_withbadwheel(ROIsR65H10, whatdoyouwant)

%%R65H10 chrimson
recording_list = {'2019_11_12_0000.mat', '2019_11_12_0001.mat', '2019_11_12_0003.mat'};



% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).leftcell = str2num(ROIsR65H10(i,1));
    ROI(i).rightcell = str2num(ROIsR65H10(i,2));

   
    
end


ROI_list = [];
rec_list = [];
cnt = 0;

if(strcmp(whatdoyouwant,'cell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).leftcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcell(j);
        end
          for j = 1:1:length(ROI(i).rightcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcell(j);
        end
      
    end
end




ROI_names =ROI;

end