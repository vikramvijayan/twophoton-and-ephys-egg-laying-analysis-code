% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS01048(ROIsSS01048, whatdoyouwant)

%function [recording] = make_small_recording_group()

% recording_list = {'2019_06_10_0000.mat','2019_06_10_0001.mat','2019_06_10_0002.mat','2019_06_10_0003.mat','2019_06_10_0005.mat',...
%     '2019_06_11_0000.mat','2019_06_11_0001.mat','2019_06_11_0002.mat','2019_06_11_0003.mat','2019_06_11_0004.mat','2019_06_11_0005.mat','2019_06_11_0006.mat',...
%     '2019_06_12_0000.mat','2019_06_12_0001.mat','2019_06_12_0002.mat','2019_06_12_0003.mat',...
%     '2019_06_13_0000.mat','2019_06_13_0001.mat','2019_06_13_0002.mat','2019_06_13_0003.mat','2019_06_13_0004.mat',...
%     '2019_06_14_0000.mat','2019_06_14_0001.mat','2019_06_14_0002.mat','2019_06_14_0003.mat',...
%     '2019_06_19_0000.mat','2019_06_19_0001.mat','2019_06_19_0002.mat','2019_06_19_0003.mat','2019_06_19_0004.mat',...
%     '2019_06_24_0000.mat','2019_06_24_0001.mat','2019_06_24_0002.mat','2019_06_24_0003.mat',...
%     '2019_06_25_0000.mat','2019_06_25_0001.mat',...
%     '2019_06_26_0000.mat','2019_06_26_0001.mat','2019_06_26_0002.mat','2019_06_26_0003.mat'};

%%SS01048
%% includes SS01048 SS53469 dual imaging
recording_list = {'2018_07_23_0000_DLC.mat','2018_07_23_0002_DLC.mat','2018_07_23_0003_DLC.mat',...
'2018_07_24_0000_DLC.mat','2018_07_24_0001_DLC.mat','2018_07_24_0002_DLC.mat','2018_07_24_0003_DLC.mat','2018_07_24_0004_DLC.mat',...
'2018_07_25_0000_DLC_new_registration.mat','2018_07_25_0001_DLC_new_registration.mat',...
'2018_09_25_0000_DLC_new_registration.mat','2018_09_25_0001_DLC_new_registration.mat','2018_09_25_0002_DLC_new_registration.mat','2018_09_25_0003_DLC_new_registration.mat','2018_09_25_0004_DLC.mat','2018_09_25_0005_DLC_new_registration.mat','2018_09_25_0006_DLC_new_registration.mat',...
'2018_09_26_0000_DLC_new_registration.mat','2018_09_26_0001_DLC_new_registration.mat','2018_09_26_0002_DLC_new_registration.mat',...
'2019_11_11_0000.mat','2019_11_11_0001.mat','2019_11_11_0002.mat','2019_11_11_0003.mat',...
'2019_11_18_0001.mat',...
'2019_11_22_0004.mat','2019_11_22_0005.mat',...
'2019_11_25_0000.mat','2019_11_25_0001.mat','2019_11_25_0002.mat','2019_11_25_0003.mat','2019_11_25_0004.mat','2019_11_25_0005.mat',...
'2019_12_17_0003.mat',...
'2019_12_18_0000.mat','2019_12_18_0001.mat','2019_12_18_0003.mat','2019_12_18_0004.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).cell = str2num(ROIsSS01048(i,1));
    ROI(i).neurite = str2num(ROIsSS01048(i,2));
    ROI(i).rightneck = str2num(ROIsSS01048(i,3));
    ROI(i).leftneck = str2num(ROIsSS01048(i,4));
    ROI(i).meanofneck = str2num(ROIsSS01048(i,6));
    ROI(i).meanofcells = str2num(ROIsSS01048(i,7));

    
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
        for j = 1:1:length(ROI(i).cell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).cell(j);
        end
      
    end
end


if(strcmp(whatdoyouwant,'leftneck'))
    for i = 1:1:length(recording_list)

        for j = 1:1:length(ROI(i).leftneck)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftneck(j);
        end
     
    end
end

if(strcmp(whatdoyouwant,'rightneck'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightneck)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightneck(j);
        end
      
     
    end
end


if(strcmp(whatdoyouwant,'neurite'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).neurite)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).neurite(j);
        end
      
      
    end
end

if(strcmp(whatdoyouwant,'meanofneck'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).meanofneck)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).meanofneck(j);
        end
      
      
    end
end

if(strcmp(whatdoyouwant,'meanofcells'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).meanofcells)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).meanofcells(j);
        end
      
      
    end
end

if(strcmp(whatdoyouwant,'anyneckside'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightneck)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightneck(j);
        end
      
        for j = 1:1:length(ROI(i).leftneck)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftneck(j);
        end   
    end
end


ROI_names =ROI;

end