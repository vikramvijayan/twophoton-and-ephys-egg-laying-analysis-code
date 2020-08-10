% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS53469_SS01048(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()

recording_list = {'2019_11_11_0000.mat',...
'2019_11_11_0001.mat',...
'2019_11_11_0002.mat',...
'2019_11_11_0003.mat',...
'2019_11_18_0001.mat',...
'2019_11_22_0004.mat',...
'2019_11_22_0005.mat',...
'2019_11_25_0000.mat',...
'2019_11_25_0001.mat',...
'2019_11_25_0002.mat',...
'2019_11_25_0003.mat',...
'2019_11_25_0004.mat',...
'2019_11_25_0005.mat',...
'2019_12_17_0003.mat',...
'2019_12_18_0000.mat',...
'2019_12_18_0001.mat',...
'2019_12_18_0003.mat',...
'2019_12_18_0004.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).leftcell = str2num(ROIsSS53469(i,1));
    ROI(i).rightcell = str2num(ROIsSS53469(i,2));
    ROI(i).leftcelldim = str2num(ROIsSS53469(i,3));
    ROI(i).rightcelldim = str2num(ROIsSS53469(i,4));
    ROI(i).dndleft = str2num(ROIsSS53469(i,5));
    ROI(i).dndright = str2num(ROIsSS53469(i,6));

    
end


ROI_list = [];
rec_list = [];
cnt = 0;

if(strcmp(whatdoyouwant,'bright_cell_ipsi'))  
    for i = 1:1:length(recording_list)
        if( (length(ROI(i).leftcell) > 0) && (length(ROI(i).dndleft) >0)) 
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcell(1);
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).dndleft(1);
        end

        if( (length(ROI(i).rightcell) > 0) && (length(ROI(i).dndright) >0)) 
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcell(1);
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).dndright(1);
        end
    end

 
end

if(strcmp(whatdoyouwant,'oviDNcell'))
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

if(strcmp(whatdoyouwant,'dndcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).dndleft)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).dndleft(j);
        end
        for j = 1:1:length(ROI(i).dndright)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).dndright(j);
        end
       
    end
end


if(strcmp(whatdoyouwant,'only_one_oviDNcell'))
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



if(strcmp(whatdoyouwant,'anyleftoviDNcell'))
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

if(strcmp(whatdoyouwant,'anyrightoviDNcell'))
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


if(strcmp(whatdoyouwant,'anybrightoviDNcell'))
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


if(strcmp(whatdoyouwant,'anydimoviDNcell'))
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


ROI_names =ROI;

end