% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS01048_chrimson_withbadwheel(ROIsSS01048, whatdoyouwant)

%function [recording] = make_small_recording_group()

%  SS01048 chrimson
recording_list = {'2019_12_11_0000.mat','2019_12_11_0001.mat','2019_12_11_0003.mat','2019_12_11_0004.mat'};


% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).cell = str2num(ROIsSS01048(i,1));
    ROI(i).neurite = str2num(ROIsSS01048(i,2));
    ROI(i).rightneck = str2num(ROIsSS01048(i,3));
    ROI(i).leftneck = str2num(ROIsSS01048(i,4));
   
    
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