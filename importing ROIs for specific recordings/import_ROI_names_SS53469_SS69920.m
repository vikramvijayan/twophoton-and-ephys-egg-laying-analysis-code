% this code is meant to help get the right ROIs for anaysis
% ROIsSS53469 is the imported table (from xlx using string)
function [rec_list, ROI_list] = import_ROI_names_SS53469_SS69920(ROIsSS53469, whatdoyouwant)

%function [recording] = make_small_recording_group()

recording_list = {'2019_12_12_0000.mat',...
'2019_12_12_0001.mat',...
'2019_12_12_0002.mat',...
'2019_12_12_0003.mat',...
'2019_12_12_0004.mat'};

% imported ROIsSS53469 as string

for i = 1:1:length(recording_list)
    ROI(i).leftcell = str2num(ROIsSS53469(i,1));
    ROI(i).rightcell = str2num(ROIsSS53469(i,2));
    ROI(i).leftcelldim = str2num(ROIsSS53469(i,3));
    ROI(i).rightcelldim = str2num(ROIsSS53469(i,4));
    ROI(i).oviENleft = str2num(ROIsSS53469(i,5));
    ROI(i).oviENleftdim = str2num(ROIsSS53469(i,6));
    ROI(i).oviENright = str2num(ROIsSS53469(i,7));
    ROI(i).oviENrightdim = str2num(ROIsSS53469(i,8));
    ROI(i).oviENneuropil = str2num(ROIsSS53469(i,9));
    
end


ROI_list = [];
rec_list = [];
cnt = 0;



if(strcmp(whatdoyouwant,'ipsi_brightcell'))  
    for i = 1:1:length(recording_list)
        if( (length(ROI(i).leftcell) > 0) && (length(ROI(i).oviENleft) >0)) 
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).leftcell(1);
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENleft(1);
        end

        if( (length(ROI(i).rightcell) > 0) && (length(ROI(i).oviENright) >0)) 
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).rightcell(1);
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENright(1);
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

if(strcmp(whatdoyouwant,'oviENcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).oviENleft)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENleft(j);
        end
        for j = 1:1:length(ROI(i).oviENright)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENright(j);
        end
        for j = 1:1:length(ROI(i).oviENleftdim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENleftdim(j);
        end
        for j = 1:1:length(ROI(i).oviENrightdim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENrightdim(j);
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

if(strcmp(whatdoyouwant,'oviENneuropil'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).oviENneuropil)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENneuropil(j);
        end
      
      
    end
end

if(strcmp(whatdoyouwant,'anybrightoviENcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).oviENright)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENright(j);
        end
      
        for j = 1:1:length(ROI(i).oviENleft)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENleft(j);
        end
     
    end
end

if(strcmp(whatdoyouwant,'anydimoviENcell'))
    for i = 1:1:length(recording_list)
        for j = 1:1:length(ROI(i).oviENrightdim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENrightdim(j);
        end
      
        for j = 1:1:length(ROI(i).oviENleftdim)
            cnt = cnt+1;
            rec_list{cnt} = char(recording_list(i));
            ROI_list(cnt) = ROI(i).oviENleftdim(j);
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