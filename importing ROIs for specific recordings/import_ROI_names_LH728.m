% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_LH728(ROIsLH728, whatdoyouwant)


    
%%LH728
recording_list = {'2019_04_29_0008_DLC.mat',...
    '2019_04_29_0009_DLC.mat',...
    '2019_04_30_0000_DLC.mat',...
    '2019_04_30_0001_DLC.mat',...
    '2019_04_30_0002_DLC.mat',...
    '2019_04_30_0003_DLC.mat'};
   



% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
 
    
    ROI(i).neurites_super = str2num(ROIsLH728(i,1));
    ROI(i).neurites_deep = str2num(ROIsLH728(i,2));
    
    
    
    %
%     if(isempty(ROI(i).leftcell))
%         ROI(i).leftcell = 0;
%     end
%     if(isempty(ROI(i).rightcell))
%         ROI(i).rightcell = 0;
%     end
%     if(isempty(ROI(i).leftcelldim))
%         ROI(i).leftcelldim = 0;
%     end
%     if(isempty(ROI(i).rightcelldim))
%         ROI(i).rightcelldim = 0;
%     end
%     if(isempty(ROI(i).leftneurite))
%         ROI(i).leftneurite = 0;
%     end
%     if(isempty(ROI(i).rightneurite))
%         ROI(i).rightneurite = 0;
%     end
%     if(isempty(ROI(i).leftneck))
%         ROI(i).leftneck = 0;
%     end
%     if(isempty(ROI(i).rightneck))
%         ROI(i).rightneck = 0;
%     end
%     if(isempty(ROI(i).combinedneck))
%         ROI(i).combinedneck = 0;
%     end
%     
end


ROI_list = [];
rec_list = [];
cnt = 0;

if(strcmp(whatdoyouwant,'anyneurite'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).neurites_super)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neurites_super(j);
        end
           for j = 1:1:length(ROI(i).neurites_deep)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neurites_deep(j);
           end
        
      
    end
end




ROI_names =ROI;

end