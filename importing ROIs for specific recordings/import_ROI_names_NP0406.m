% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_NP0406(ROIsNP0406, whatdoyouwant)


recording_list = {'2017_02_21_0002N_DLC.mat',...
    '2017_08_02_0003N_DLC.mat',...
    '2017_08_02_0002N_DLC.mat',...
    '2017_08_22_0003N_DLC.mat',...
    '2017_08_22_0003N_A_DLC.mat',...
    '2017_08_22_0003N_B_DLC.mat',...
    '2017_08_23_0000N_DLC.mat',...
    '2017_11_03_0000N_DLC.mat',...
    '2017_12_05_0000N_DLC.mat',...
    '2017_12_05_0001N_DLC.mat',...
    '2017_12_05_0001N_A_DLC.mat',...
    '2017_12_06_0002N_DLC.mat',...
    '2017_12_06_0002N_A_DLC.mat',...
    '2017_12_06_0002N_B_DLC.mat',...
    '2017_12_28_0000N_DLC.mat',...
    '2017_12_28_0003N_DLC.mat',...
    '2017_12_28_0004N_DLC.mat',...
    '2017_12_31_0000N_DLC.mat',...
    '2017_12_31_0001N_DLC.mat',...
    '2018_01_02_0000N_DLC.mat',...
    '2018_01_10_0000_DLC.mat',...
    '2018_01_10_0001_DLC.mat',...
    '2018_01_15_0000_DLC.mat',...
    '2018_01_15_0001_DLC.mat',...
    '2018_01_15_0002_DLC.mat',...
    '2018_01_15_0003_DLC.mat',...
    '2018_01_16_0000_DLC.mat',...
    '2018_07_25_0002_DLC.mat',...
    '2018_07_25_0003_DLC.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).neurite = str2num(ROIsNP0406(i,1));
    ROI(i).single_cell_close = str2num(ROIsNP0406(i,2));
    
    ROI(i).cell_left = str2num(ROIsNP0406(i,3));
    ROI(i).cell_right = str2num(ROIsNP0406(i,4));
    ROI(i).neck_left = str2num(ROIsNP0406(i,5));
    ROI(i).neck_right = str2num(ROIsNP0406(i,6));
    
    ROI(i).dimcell = str2num(ROIsNP0406(i,7));
    
    
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