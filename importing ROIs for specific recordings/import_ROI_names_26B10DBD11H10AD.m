% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_26B10DBD11H10AD(ROIs26B10DBD11H10AD, whatdoyouwant)

%%NP0406 split  26B10DBD 11H10AD
recording_list = { '2018_09_26_0003_DLC.mat','2018_09_26_0004_DLC.mat',...
    '2018_09_27_0000_DLC.mat', '2018_09_27_0001_DLC.mat', '2018_09_27_0002_DLC.mat', '2018_09_27_0003_DLC.mat', '2018_09_27_0004_DLC.mat',...
    '2018_10_10_0000_DLC.mat', '2018_10_10_0001_DLC.mat', '2018_10_10_0002_DLC.mat', '2018_10_10_0003_DLC.mat', '2018_10_10_0004_DLC.mat', '2018_10_10_0005_DLC.mat',...
    '2018_10_11_0000_DLC.mat','2018_10_11_0001_DLC.mat','2018_10_11_0002_DLC.mat','2018_10_11_0003_DLC.mat','2018_10_11_0004_DLC.mat','2018_10_11_0005_DLC.mat','2018_10_11_0006_DLC.mat',...
    '2018_10_12_0000_DLC.mat', '2018_10_12_0001_DLC.mat', '2018_10_12_0002_DLC.mat', '2018_10_12_0003_DLC.mat', '2018_10_12_0004_DLC.mat'};


% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).neurite = str2num(ROIs26B10DBD11H10AD(i,1));
    ROI(i).single_cell_close = str2num(ROIs26B10DBD11H10AD(i,2));
    
    ROI(i).cell_left = str2num(ROIs26B10DBD11H10AD(i,3));
    ROI(i).cell_right = str2num(ROIs26B10DBD11H10AD(i,4));
    ROI(i).neck_left = str2num(ROIs26B10DBD11H10AD(i,5));
    ROI(i).neck_right = str2num(ROIs26B10DBD11H10AD(i,6));
    
    ROI(i).dimcell = str2num(ROIs26B10DBD11H10AD(i,7));
    
    
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

if(strcmp(whatdoyouwant,'cell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).cell_left)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell_left(j);
        end
           for j = 1:1:length(ROI(i).cell_right)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell_right(j);
        end
      
    end
end


if(strcmp(whatdoyouwant,'dimcell'))
    for i = 1:1:length(recording_list)

        for j = 1:1:length(ROI(i).dimcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).dimcell(j);
        end
     
    end
end

if(strcmp(whatdoyouwant,'singlecellclose'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).single_cell_close)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).single_cell_close(j);
        end
      
     
    end
end


if(strcmp(whatdoyouwant,'anyneurite'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).neurite)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neurite(j);
        end
    end
end

if(strcmp(whatdoyouwant,'anyneckside'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).neck_right)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neck_right(j);
        end
      
        for j = 1:1:length(ROI(i).neck_left)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neck_left(j);
        end   
    end
end


ROI_names =ROI;

end