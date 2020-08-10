% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS53469_hardsoft(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()

recording_list = {'2019_11_14_0001.mat','2019_11_14_0002.mat','2019_11_14_0003.mat',...
'2019_12_13_0000.mat',...
'2019_12_13_0001.mat',...
'2019_12_13_0002.mat',...
'2019_12_13_0003.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).leftcell = str2num(ROIsSS53469(i,1));
    ROI(i).rightcell = str2num(ROIsSS53469(i,2));
    ROI(i).leftcelldim = str2num(ROIsSS53469(i,3));
    ROI(i).rightcelldim = str2num(ROIsSS53469(i,4));
    ROI(i).leftneurite = str2num(ROIsSS53469(i,5));
    ROI(i).rightneurite = str2num(ROIsSS53469(i,6));
    ROI(i).leftneck = str2num(ROIsSS53469(i,7));
    ROI(i).rightneck = str2num(ROIsSS53469(i,8));
    ROI(i).combinedneck = str2num(ROIsSS53469(i,9));
    
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
        for j = 1:1:length(ROI(i).leftcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcelldim(j);
        end
        for j = 1:1:length(ROI(i).rightcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcelldim(j);
        end
    end
end


if(strcmp(whatdoyouwant,'only_one_cell'))
    for i = 1:1:length(recording_list)
        cell_taken = 0;
        for j = 1:1:length(ROI(i).leftcell)
            if(cell_taken == 0)
                cnt = cnt+1;
                rec_list{cnt} = char(recording_list(i));
                ROI_list(cnt) = ROI(i).leftcell(j);
                cell_taken=1;
            end
        end
        for j = 1:1:length(ROI(i).rightcell)
            if(cell_taken == 0)
                cnt = cnt+1;
                rec_list{cnt} = char(recording_list(i));
                ROI_list(cnt) = ROI(i).rightcell(j);
                cell_taken=1;
            end
        end
        for j = 1:1:length(ROI(i).leftcelldim)
            if(cell_taken == 0)
                cnt = cnt+1;
                rec_list{cnt} = char(recording_list(i));
                ROI_list(cnt) = ROI(i).leftcelldim(j);
                cell_taken=1;
            end
        end
        for j = 1:1:length(ROI(i).rightcelldim)
            if(cell_taken == 0)
                cnt = cnt+1;
                rec_list{cnt} = char(recording_list(i));
                ROI_list(cnt) = ROI(i).rightcelldim(j);
                cell_taken=1;
            end
        end
    end
end



if(strcmp(whatdoyouwant,'anyleftcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).leftcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcell(j);
        end
      
        for j = 1:1:length(ROI(i).leftcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcelldim(j);
        end
     
    end
end

if(strcmp(whatdoyouwant,'anyrightcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcell(j);
        end
      
        for j = 1:1:length(ROI(i).rightcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcelldim(j);
        end
     
    end
end


if(strcmp(whatdoyouwant,'anybrightcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcell(j);
        end
      
        for j = 1:1:length(ROI(i).leftcell)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcell(j);
        end
     
    end
end

if(strcmp(whatdoyouwant,'anydimcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcelldim(j);
        end
      
        for j = 1:1:length(ROI(i).leftcelldim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcelldim(j);
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

if(strcmp(whatdoyouwant,'anyneurite'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).rightneurite)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightneurite(j);
        end
      
        for j = 1:1:length(ROI(i).leftneurite)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftneurite(j);
        end   
    end
end

ROI_names =ROI;

end