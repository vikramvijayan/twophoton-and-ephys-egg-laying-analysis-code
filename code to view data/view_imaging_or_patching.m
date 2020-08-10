% put tseries  -1 for patch (second input)
function varargout = view_imaging_or_patching(varargin)
% VIEW_IMAGING_OR_PATCHING MATLAB code for view_imaging_or_patching.fig
%      VIEW_IMAGING_OR_PATCHING, by itself, creates a new VIEW_IMAGING_OR_PATCHING or raises the existing
%      singleton*.
%
%      H = VIEW_IMAGING_OR_PATCHING returns the handle to a new VIEW_IMAGING_OR_PATCHING or the handle to
%      the existing singleton*.
%
%      VIEW_IMAGING_OR_PATCHING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_IMAGING_OR_PATCHING.M with the given input arguments.
%
%      VIEW_IMAGING_OR_PATCHING('Property','Value',...) creates a new VIEW_IMAGING_OR_PATCHING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_imaging_or_patching_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_imaging_or_patching_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_imaging_or_patching

% Last Modified by GUIDE v2.5 08-Oct-2019 18:40:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @view_imaging_or_patching_OpeningFcn, ...
    'gui_OutputFcn',  @view_imaging_or_patching_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
warning('off');

% End initialization code - DO NOT EDIT
end

% --- Executes just before view_imaging_or_patching is made visible.
function view_imaging_or_patching_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_imaging_or_patching (see VARARGIN)


warning('off');
%matlab.video.read.UseHardwareAcceleration('off')
matlab.video.read.UseHardwareAcceleration('on')

recording = varargin{1};
tseries_for_plot = varargin{2};

if(tseries_for_plot ==-1)
    handles.patch = 1;
else
    handles.patch = 0;
end
%varargin 3 is useless
%if(varargin{3} == 1)
handles.filename = recording.movie1.filename;
%else
handles.filename2 = recording.movie2.filename;
%end


if(length(varargin) > 3)
    handles.time_st = varargin{4};
    handles.time_en = varargin{5};
else
    if(handles.patch == 0)
        handles.time_st = recording.tseries(tseries_for_plot).Time_s(1);
        handles.time_en = recording.tseries(tseries_for_plot).Time_s(end);
    else
        handles.time_st = recording.time_to_use(1);
        handles.time_en = recording.time_to_use(2);
    end
end

if(handles.patch == 0)
    
    InputStack = recording.tseries(tseries_for_plot).tsStack;
end
set(handles.fps, 'String', '25' );
handles.stop_play = 1;

if(strcmp(handles.filename((end-3):(end)),'.fmf'))
    % extract header from the appropriate FMF
    [header_size, version, f_height, f_width, bytes_per_chunk, max_n_frames, data_format] = fmf_read_header(handles.filename);
    handles.header_size = header_size;
    handles.version = version;
    handles.f_height = f_height;
    handles.f_width = f_width;
    handles.bytes_per_chunk = bytes_per_chunk;
    handles.max_n_frames = max_n_frames;
    handles.data_format = data_format;
    handles.fp = fopen( handles.filename, 'r' );
    %fseek(handles.fp,handles.header_size,'bof');
    %[data, stamp] = fmf_read_frame( handles.fp, handles.f_height, handles.f_width, handles.bytes_per_chunk, handles.data_format );
    handles.movie_length = max_n_frames;
    handles.isfmf = 1;
end

if(strcmp(handles.filename((end-3):(end)),'ufmf'))
    header = ufmf_read_header(handles.filename);
    handles.header = header;
    handles.movie_length = header.nframes;
    handles.max_n_frames = header.nframes;
    handles.f_height = header.max_width;
    handles.f_width = header.max_height;
    handles.isfmf = 2;
end

if(strcmp(handles.filename((end-3):(end)),'.avi') || strcmp(handles.filename((end-3):(end)),'.mp4'))
    handles.vobj = VideoReader(handles.filename);
    handles.max_n_frames = handles.vobj.NumberOfFrames;
    handles.movie_length = handles.vobj.NumberOfFrames;
    handles.f_height = handles.vobj.Height;
    handles.f_width = handles.vobj.Width;
    handles.isfmf = 0;
    delete(handles.vobj);
    handles.vobj = VideoReader(handles.filename);
    
    handles.vobj2 = VideoReader(handles.filename2);
    handles.max_n_frames2 = handles.vobj2.NumberOfFrames;
    handles.movie_length2 = handles.vobj2.NumberOfFrames;
    handles.f_height2 = handles.vobj2.Height;
    handles.f_width2 = handles.vobj2.Width;
    delete(handles.vobj2);
    handles.vobj2 = VideoReader(handles.filename2);
    
end

handles.recording = recording;
handles.tseries_for_plot = tseries_for_plot;

if(handles.patch ==0)
    handles.tsStack = InputStack;
    size_mat = size(InputStack);
    handles.size_mat = size_mat;
    handles.total_ROI = 0;
    
end
if(handles.patch ==1)
    handles.size_mat = [0,0,0,0,2];
    size_mat = [0,0,0,0,2];
end


% if there is more than 1 z slice
if(length(size_mat) == 5)
    handles.ROI_masks = zeros(size_mat(1), size_mat(2), size_mat(5),20);
    
    set(handles.z_slider,'Max',handles.size_mat(5),'Min',1);
    set(handles.z_slider,'SliderStep',[1/(size_mat(5)-1) 1/(size_mat(5)-1)]);
    set(handles.z_slider, 'Value', 1);
    
    set(handles.t_slider,'Max',handles.movie_length,'Min',1);
    set(handles.t_slider,'SliderStep',[1/(handles.movie_length-1) 100/(handles.movie_length-1)]);
    set(handles.t_slider, 'Value', 1);
else
    handles.ROI_masks = zeros(size_mat(1), size_mat(2), 1,20);
    
    set(handles.z_slider,'Max',1,'Min',1);
    set(handles.z_slider,'SliderStep',[0 0]);
    set(handles.z_slider, 'Value', 1);
    
    set(handles.t_slider,'Max',handles.movie_length,'Min',1);
    set(handles.t_slider,'SliderStep',[10/(handles.movie_length-1) 100/(handles.movie_length-1)]);
    set(handles.t_slider, 'Value', 1);
end

set(handles.maxcolor,'Max',2^12,'Min',0);
set(handles.maxcolor,'SliderStep',[1/(2^12-1) 10/(2^12-1)]);
set(handles.maxcolor, 'Value', 2^12);

set(handles.mincolor,'Max',2^12,'Min',0);
set(handles.mincolor,'SliderStep',[1/(2^12-1) 10/(2^12-1)]);
set(handles.mincolor, 'Value', 0);


set(handles.slidermax,'Max',2^8,'Min',0);
set(handles.slidermax,'SliderStep',[1/(2^8-1) 10/(2^8-1)]);
set(handles.slidermax, 'Value', 2^8);

set(handles.slidermin,'Max',2^8,'Min',0);
set(handles.slidermin,'SliderStep',[1/(2^8-1) 10/(2^8-1)]);
set(handles.slidermin, 'Value', 0);


%set(handles.sucrose_roi_thresh,'Max',max(recording.movie1.roisucrose),'Min',min(recording.movie1.roisucrose));
%set(handles.sucrose_roi_thresh,'SliderStep',[1/(max(recording.movie1.roisucrose)-1) 10/(max(recording.movie1.roisucrose)-1)]);
%set(handles.sucrose_roi_thresh, 'Value', median(recording.movie1.roisucrose));

set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
set(handles.current_z_display, 'String', num2str(ceil(handles.z_slider.Value)));

set(handles.cmap_select, 'String', {'gray','parula','jet','hsv','hot','cool'});
handles.cmap_select.Value = 1;
handles.cmap_options  =  {'gray','parula','jet','hsv','hot','cool'};



handles.first_run = 1;

% was having some problems with scroll wheel so disabled for the moment
% set(gcf, 'WindowScrollWheelFcn', {@wheel,handles});
% handles.sliderListener = addlistener(handles.t_slider,'ContinuousValueChange', @(hObject, event) t_slider_Callback(hObject, eventdata, handles));

% Choose default command line output for view_imaging_or_patching
handles.output = hObject;

% plot the imaging data
if(handles.patch ==0)
    tseries = handles.recording.tseries(handles.tseries_for_plot);
    
    axes(handles.axes3);  hold on;
    %imagesc(tseries.Time_s,1:tseries.ROI_number,(tseries.df_over_f)./nanmean(tseries.df_over_f,2)   ); colormap(gca,'jet');
    interpsignal = interp1(tseries.Time_s,tseries.df_over_f(1,:)./nanmean(tseries.df_over_f(1,:)), tseries.Time_s(1):.1:tseries.Time_s(end),'previous');
    interpsignal_smooth = smoothdata(interpsignal,'movmean',50,'includenan');
    
    
    
    
    
    plot(tseries.Time_s(1):.1:tseries.Time_s(end),interpsignal_smooth,'k');
    
    
    
    %     imagesc(tseries.Time_s(1):.1:tseries.Time_s(end),1:1,interpsignal);
    %     set(gca,'ytick',[])
    %     m=1000;
    %     cm_magma=magma(m);
    %     cm_inferno=inferno(m);
    %     cm_plasma=plasma(m);
    %     cm_viridis=viridis(m);
    %     colormap(gca,cm_plasma);
    %     %colormap(gca,'jet');
    %     caxis([-1,4]);
    
    %imagesc(tseries.Time_s,1:tseries.ROI_number,(tseries.df_over_f) ); colormap(gca,'jet'); caxis([.1,1.05]);
    
    %ylabel('ROIs','fontsize', 8);
    set(gca,'xlim',[handles.time_st, handles.time_en]);
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.004 0.004]);
    set(gca,'FontSize',8)
    %     set(gca,'ylim',[0.5 tseries.ROI_number+.5]);
    set(gca,'ylim',[-.1 4.1]);
    
    %     set(gca,'ytick',[])
    box on;
end

% plot patching
if(handles.patch ==1)
    %tseries = handles.recording.tseries(handles.tseries_for_plot);
    
    axes(handles.axes3);  hold on;
    plot(handles.recording.abf.Time_s(1:100:end), 10000.*handles.recording.abf.CH1_patch_spikes_conv_rect(1:100:end))
    ylabel('ROIs','fontsize', 8);
    set(gca,'ylim',[0 (max(10000.*handles.recording.abf.CH1_patch_spikes_conv_rect)+5)]);
    set(gca,'xlim',[handles.time_st, handles.time_en]);
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.004 0.004]);
    set(gca,'FontSize',8)
    box on;
end


% plot the trajectory or wheel data
axes(handles.axes5);  hold on;
% just added
%yyaxis right;
% plot(recording.movie1.time_stamps,recording.movie1.abd_length./600,'-k');
% plot(recording.movie1.time_stamps,recording.movie1.abd_angle,'-m');
%set(gca,'ylim',[.3 1.1]);

%set(gca,'ytick',[]);

%yyaxis left;
%yyaxis left;

[a b] = find(recording.movie1.sucrose == 200);
%scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'g');
scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'MarkerEdgeColor',[33, 113, 181]./255,'MarkerFaceColor',[33, 113, 181]./255);


[a b] = find(recording.movie1.sucrose == 0);
%scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'r');
scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'MarkerEdgeColor',[158, 202, 225]./255,'MarkerFaceColor',[158, 202, 225]./255);


[a b] = find(recording.movie1.sucrose == 500);
%scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'b');
scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'MarkerEdgeColor',[8, 48, 107]./255,'MarkerFaceColor',[8, 48, 107]./255);

[a b] = find(recording.movie1.sucrose == 1);
%scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'g');
scatter(handles.recording.movie1.time_stamps(b),recording.movie1.filtered_wheel(b) ,.5,'MarkerEdgeColor',[8, 48, 107]./255,'MarkerFaceColor',[8, 48, 107]./255);


scatter(handles.recording.movie1.time_stamps(handles.recording.movie1.eggs), recording.movie1.filtered_wheel(handles.recording.movie1.eggs),25,'r','filled');

line([handles.recording.movie1.time_stamps(handles.recording.movie1.eggs), handles.recording.movie1.time_stamps(handles.recording.movie1.eggs)], [-2 6],'color','r');
scatter([handles.recording.movie1.time_stamps(handles.recording.movie1.eggs)], ones(1,length(handles.recording.movie1.eggs)).*-1,25,'r');

%behavior = [16075,17948,18348,20082,20270,20757,20761,20944;44292,44980,45800,46730,47075,47901,47982,48156;54738,55232,55736,56670,56760,57714,57773,57900;72718,73151,73655,74164,74168,75034,75560,75670;64070,64685,65233,66042,66054,66470,66963,67193;96344,96857,97318,98103,98911,99425,99769,99945;102472,103100,103613,104478,104605,105598,105612,105879;108118,108464,109095,109886,110134,110619,110619,110893;112893,113110,113638,114327,114693,116085,116117,116299;118694,118884,119896,120790,120877,122120,122120,122276;139915,140472,141146,141644,142421,143266,143266,143390];


behavior = [29242,29560,30086,30770,30770,31805,32093,32341;44415,44956,45383,45800,45800,47360,47455,47606;76005,79248,79917,80413,80431,82633,82682,82727;159234,159248,160253,160979,161591,162511,162622,162709];
    


[numlines, ~] = size(behavior);
for i = 1:1:numlines
    line([handles.recording.movie1.time_stamps(behavior(i,1)), handles.recording.movie1.time_stamps(behavior(i,1))], [-2 16],'color',[216,191,216]./255);
    scatter([handles.recording.movie1.time_stamps(behavior(i,1))], recording.movie1.filtered_wheel(behavior(i,1)),25,[216,191,216]./255,'filled');
    
    
    line([handles.recording.movie1.time_stamps(behavior(i,6)), handles.recording.movie1.time_stamps(behavior(i,6))],[-2 16], 'color',[128,0,128]./255);
    scatter([handles.recording.movie1.time_stamps(behavior(i,6))], recording.movie1.filtered_wheel(behavior(i,6)),25,[128,0,128]./255,'filled');
    
    line([handles.recording.movie1.time_stamps(behavior(i,7)), handles.recording.movie1.time_stamps(behavior(i,7))],[-2 16], 'color','r');
    scatter([handles.recording.movie1.time_stamps(behavior(i,7))], [-1],25,'r','filled');
end


set(gca,'TickDir','out');
set(gca,'xlim',[handles.recording.movie1.time_stamps(1), handles.recording.movie1.time_stamps(end)]);
set(gca,'ylim',[0, 2*pi]);
set(gca,'FontSize',8);
set(gca,'TickLength',[0.004 0.004]);
set(gca,'xtick',[])

%yyaxis right; hold on;
%plot(handles.recording.abf.Time_s, handles.recording.abf.PWMlaser_uWpermm2,'m');

box on;

linkaxes([handles.axes5,handles.axes3],'x');

% save the changes to handles
guidata(hObject, handles);

% UIWAIT makes view_imaging_or_patching wait for user response (see UIRESUME)
end

% --- Outputs from this function are returned to the command line.
function varargout = view_imaging_or_patching_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
warning('on');
%delete(handles.figure1);
end



% --- Executes on slider movement.
function z_slider_Callback(hObject, eventdata, handles)
% hObject    handle to z_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
set(handles.current_z_display, 'String', num2str(ceil(handles.z_slider.Value)));
guidata(hObject, handles);

populate_graphs(hObject);



end

% --- Executes during object creation, after setting all properties.
function z_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function t_slider_Callback(hObject, eventdata, handles)
% hObject    handle to t_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
set(handles.current_z_display, 'String', num2str(ceil(handles.z_slider.Value)));
guidata(hObject, handles);

populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function t_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function populate_graphs(hObject)

handles = guidata(hObject);

plot_data_line(hObject);

if(handles.patch == 0)
    plot_ROI(hObject);
end

if(handles.isfmf == 0)
    plot_movie_avi(hObject);
end

if(handles.isfmf == 1)
    plot_movie_fmf(hObject);
end

if(handles.isfmf == 2)
    plot_movie_ufmf(hObject);
end

drawnow;
end

function plot_ROI(hObject)
handles = guidata(hObject);

map = colormap(gca,char(handles.cmap_options(handles.cmap_select.Value)));

% find the closest 2p image to the behavior
[~, closest_image] = min(abs(handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))-handles.recording.tseries(handles.tseries_for_plot).Time_s));
set(handles.abf_time, 'String', num2str(handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))));


% if there are z
if(length(handles.size_mat) == 5)
    
    if(~handles.max_t.Value && ~handles.max_z.Value)
        curr_image = handles.tsStack(:,:,1,closest_image, round(handles.z_slider.Value));
    end
    
    if(handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,1,:,round(handles.z_slider.Value)),[],4);
    end
    
    if(handles.mean_z.Value)
        curr_image = mean(handles.tsStack(:,:,1,closest_image,:),5);
    end
    
    if(handles.max_z.Value)
        curr_image = max(handles.tsStack(:,:,1,closest_image,:),[],5);
    end
    
    if(handles.max_z.Value & handles.max_t.Value)
        curr_image = max(max(handles.tsStack(:,:,1,:,:),[],5),[],4);
    end
    
else
    if(~handles.max_t.Value && ~handles.max_z.Value)
        curr_image = handles.tsStack(:,:,closest_image);
    end
    
    if(handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,:),[],3);
    end
    
    if(handles.max_z.Value)
        curr_image = handles.tsStack(:,:,closest_image);
    end
    
    if(handles.max_z.Value & handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,:),[],3);
    end
end



if(handles.first_run)
    axes(handles.axes1); cla('reset'); hold on;
    
    m=1000;
    cm_magma=magma(m);
    cm_inferno=inferno(m);
    cm_plasma=plasma(m);
    cm_viridis=viridis(m);
    
    %handles.ROI_show = imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)],'Colormap',cm_plasma);
    handles.ROI_show = imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)],'Colormap',jet(1000));
    
else
    axes(handles.axes1);  cla('reset'); hold on;
    delete(handles.ROI_show);
    m=1000;
    cm_magma=magma(m);
    cm_inferno=inferno(m);
    cm_plasma=plasma(m);
    cm_viridis=viridis(m);
    % handles.ROI_show = imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)],'Colormap',cm_plasma);
    handles.ROI_show = imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)],'Colormap',jet(1000));
    
    %set(handles.ROI_show ,'cdata',curr_image);
end

%colormap(map);
%colormap(gca,cm_plasma);
colormap(gca,'jet');

% put the ROI as outline
if(handles.plotselectedrois.Value == 1)
    if(handles.total_ROI == 0)
        total_to_use = handles.recording.tseries(handles.tseries_for_plot).ROI_number;
        rois_to_use = handles.recording.tseries(handles.tseries_for_plot).ROI;
    else
        total_to_use = handles.total_ROI;
        rois_to_use = handles.ROI_masks;
    end
    
    
    if(handles.max_z.Value == 0)
        
        for i = 1:1:total_to_use
            b = bwboundaries(rois_to_use(:,:,ceil(handles.z_slider.Value),i));
            if(~isempty(b))
                tt = text(b{1,1}(1,2),b{1,1}(1,1),num2str(i));
                tt.Color = 'green';
                for k = 1:numel(b)
                    plot(b{k}(:,2), b{k}(:,1), 'r', 'Linewidth', .5);
                end
            end
        end
    else
        for j =1:1:handles.size_mat(5)
            for i = 1:1:total_to_use
                b = bwboundaries(rois_to_use(:,:,j,i));
                tt = text(b{1,1}(1,2),b{1,1}(1,1),num2str(i));
                tt.Color = 'green';
                for k = 1:numel(b)
                    plot(b{k}(:,2), b{k}(:,1), 'r', 'Linewidth', .5);
                end
            end
        end
    end
end

guidata(hObject, handles);

end


function plot_movie_fmf(hObject)
handles = guidata(hObject);


fseek(handles.fp,handles.header_size+(ceil(handles.t_slider.Value)-1)*handles.bytes_per_chunk,'bof');
[data, stamp] = fmf_read_frame( handles.fp, handles.f_height, handles.f_width, handles.bytes_per_chunk, handles.data_format );

% need to change this to the stamp correcte for abf start
% set(handles.moviestamp, 'String', num2str(stamp-handles.first_stamp));

vidFrame = double(data);

% Commented from reading avi
% handles.recording.movie1.videoobj.CurrentTime = ceil(handles.t_slider.Value)./handles.recording.movie1.videoobj.FrameRate - 1/handles.recording.movie1.videoobj.FrameRate;
% vidFrame = readFrame(handles.recording.movie1.videoobj);

if(handles.first_run)
    axes(handles.axes2); hold on;
    handles.movie_show = imagesc(flipud(vidFrame));
    axis tight;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    handles.first_run = 0;
    map = colormap(gca,char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(gca, map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
else
    axes(handles.axes2); hold on;
    set(handles.movie_show ,'cdata',flipud(vidFrame));
    map = colormap(gca,char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(gca,map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end

guidata(hObject, handles);
end


function plot_movie_ufmf(hObject)
handles = guidata(hObject);

[im,header,timestamp,bb,mu] = ufmf_read_frame(handles.header,ceil(handles.t_slider.Value));


% need to change this to the stamp correcte for abf start
% set(handles.moviestamp, 'String', num2str(stamp-handles.first_stamp));

vidFrame = double(im);

% Commented from reading avi
% handles.recording.movie1.videoobj.CurrentTime = ceil(handles.t_slider.Value)./handles.recording.movie1.videoobj.FrameRate - 1/handles.recording.movie1.videoobj.FrameRate;
% vidFrame = readFrame(handles.recording.movie1.videoobj);

if(handles.first_run)
    axes(handles.axes2); hold on;
    handles.movie_show = imagesc(flipud(vidFrame));
    axis tight;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    handles.first_run = 0;
    map = colormap(gca,char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(gca,map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
else
    axes(handles.axes2); hold on;
    set(handles.movie_show ,'cdata',flipud(vidFrame));
    map = colormap(gca,char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(gca,map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end

guidata(hObject, handles);
end




function plot_movie_avi(hObject)
handles = guidata(hObject);


handles.vobj.CurrentTime = ceil(handles.t_slider.Value)./handles.vobj.FrameRate - 1/handles.vobj.FrameRate;
vidFrame_m1 = readFrame(handles.vobj);

vidFrame = double(vidFrame_m1(:,:,1));


handles.vobj2.CurrentTime = (handles.movie_length2 - (handles.movie_length - ceil(handles.t_slider.Value)))./handles.vobj2.FrameRate - 1/handles.vobj2.FrameRate;
vidFrame_m2 = readFrame(handles.vobj2);

vidFrame2 = double(vidFrame_m2(:,:,1));
clrz = 'rgbmcy';

clrz = [102, 0, 204; 153 0 255; 204 0 255; 255 0 255; 255 0 102; 113 196 120]./255;
if(handles.first_run)
    axes(handles.axes2); hold on;
    %handles.movie_show = imagesc(flipud(vidFrame));
    handles.movie_show = imshow(vidFrame,[handles.slidermin.Value, handles.slidermax.Value],'Colormap',gray(1000));
    axis tight;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    cnt = 0;
    for i =1:3:16
        cnt = cnt+1;
        index_of_frame = round((handles.vobj.CurrentTime+1/handles.vobj.FrameRate).*handles.vobj.FrameRate);
        if(handles.recording.movie1.DLC(i+2,index_of_frame) > .2)
            scatter(handles.recording.movie1.DLC(i,index_of_frame),handles.recording.movie1.DLC(i+1,index_of_frame),25,clrz(cnt,:),'filled');
        end
    end
    
    axes(handles.axes6); hold on;
    %handles.movie_show = imagesc(flipud(vidFrame));
    handles.movie_show2 = imshow(vidFrame2,[handles.slidermin.Value, handles.slidermax.Value],'Colormap',gray(1000));
    axis tight;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    cnt = 0;
    
    for i =1
        cnt = cnt+1;
        index_of_frame = round((handles.vobj2.CurrentTime+1/handles.vobj2.FrameRate).*handles.vobj2.FrameRate);
        if(handles.recording.movie2.DLC(i+2,index_of_frame) > .2)
            scatter(handles.recording.movie2.DLC(i,index_of_frame),handles.recording.movie2.DLC(i+1,index_of_frame),25,'b','filled');
        end
    end
    
    handles.first_run = 0;
    %map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    %colormap(map);
    %set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
    
else
    axes(handles.axes2); hold on; %cla('reset'); hold on;
    axesHandlesToChildObjects = findobj(handles.axes2, 'Type', 'image');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    
    axesHandlesToChildObjects = findobj(handles.axes2, 'Type', 'scatter');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    
    handles.movie_show = imshow(vidFrame,[handles.slidermin.Value, handles.slidermax.Value],'Colormap',gray(1000));
    cnt = 0;
    
    for i =1:3:16
        cnt = cnt+1;
        index_of_frame = round((handles.vobj.CurrentTime+1/handles.vobj.FrameRate).*handles.vobj.FrameRate);
        if(handles.recording.movie1.DLC(i+2,index_of_frame) > .2)
            scatter(handles.recording.movie1.DLC(i,index_of_frame),handles.recording.movie1.DLC(i+1,index_of_frame),25,clrz(cnt,:),'filled');
        end
    end
    
    axes(handles.axes6); hold on; %cla('reset'); hold on;
    axesHandlesToChildObjects = findobj(handles.axes6, 'Type', 'image');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    
    axesHandlesToChildObjects = findobj(handles.axes6, 'Type', 'scatter');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    handles.movie_show2 = imshow(vidFrame2,[handles.slidermin.Value, handles.slidermax.Value],'Colormap',gray(1000));
    cnt = 0;
    
    for i =1
        cnt=cnt+1;
        index_of_frame = round((handles.vobj2.CurrentTime+1/handles.vobj2.FrameRate).*handles.vobj2.FrameRate);
        if(handles.recording.movie2.DLC(i+2,index_of_frame) > .2)
            scatter(handles.recording.movie2.DLC(i,index_of_frame),handles.recording.movie2.DLC(i+1,index_of_frame),10,'b', 'filled');
        end
    end
    %set(handles.movie_show ,'cdata',flipud(vidFrame));
    %map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    %colormap(map);
    %set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end


guidata(hObject, handles);
end



% function wheel(hObject, callbackdata, handles)
% handles = guidata(hObject);
%
% if callbackdata.VerticalScrollCount > 0
%     handles.t_slider.Value = min(handles.movie_length,handles.t_slider.Value + 10);
% elseif callbackdata.VerticalScrollCount < 0
%     handles.t_slider.Value =max(1,handles.t_slider.Value - 10);
% end
% set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
%
% guidata(hObject,handles);
%
% populate_graphs(hObject);
% end


function plot_data_line(hObject)
handles = guidata(hObject);

if(handles.patch ==0)
    [~, closest_image] = min(abs(handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))-handles.recording.tseries(handles.tseries_for_plot).Time_s));
end
if(handles.patch==1)
    closest_image = handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))*10000;
end

if(~handles.first_run)
    axes(handles.axes3); hold on;
    children = get( gca, 'children');
    delete(children(1));
    
    axes(handles.axes5); hold on;
    children = get( gca, 'children');
    delete(children(1));
end

axes(handles.axes3);
if(handles.patch ==0)
    line([handles.recording.tseries(handles.tseries_for_plot).Time_s(closest_image),handles.recording.tseries(handles.tseries_for_plot).Time_s(closest_image)],[-10,10],'Color','m');
end
if(handles.patch==1)
    line([handles.recording.abf.Time_s(floor(closest_image)),handles.recording.abf.Time_s(floor(closest_image))],[-10,10],'Color','m');
end



axes(handles.axes5);
line([handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value)),handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))], [-10,10],'Color','m');


guidata(hObject, handles);

end


% --- Executes on button press in show_roi.
function show_roi_Callback(hObject, eventdata, handles)
% hObject    handle to show_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
populate_graphs(hObject);
% Hint: get(hObject,'Value') returns toggle state of show_roi
end


% --- Executes on button press in max_t.
function max_t_Callback(hObject, eventdata, handles)
% hObject    handle to max_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
populate_graphs(hObject);

% Hint: get(hObject,'Value') returns toggle state of max_t
end

% --- Executes on button press in max_z.
function max_z_Callback(hObject, eventdata, handles)
% hObject    handle to max_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure populate_graphs(hObject);
% with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of max_z
end


% --- Executes on selection change in cmap_select.
function cmap_select_Callback(hObject, eventdata, handles)
% hObject    handle to cmap_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap_select
populate_graphs(hObject);
end

% --- Executes during object creation, after setting all properties.
function cmap_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmap_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on slider movement.
function mincolor_Callback(hObject, eventdata, handles)
% hObject    handle to mincolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function mincolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mincolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function maxcolor_Callback(hObject, eventdata, handles)
% hObject    handle to maxcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function maxcolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on button press in zoomin.
function zoomin_Callback(hObject, eventdata, handles)
% hObject    handle to zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2); hold on;
zoom on;


end


% --- Executes on button press in zoomoff.
function zoomoff_Callback(hObject, eventdata, handles)
% hObject    handle to zoomoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2); hold on;
zoom off;
end


% --- Executes on button press in markboundary.
function markboundary_Callback(hObject, eventdata, handles)
% hObject    handle to markboundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ctime = handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value));
[a b] = min(abs(handles.recording.abf.Time_s-ctime));
set(handles.boundary, 'String', num2str(handles.recording.abf.Wheel(b)));



if(~handles.first_run)
    
    axes(handles.axes5); hold on;
    children = get( gca, 'children');
    delete(children(1));
end

axes(handles.axes5);

otherbound = mod(handles.recording.abf.Wheel(b) + 1.5,4);
if(otherbound < handles.recording.abf.Wheel(b))
    otherbound = otherbound + 1;
end

line([handles.recording.movie1.time_stamps(1),handles.recording.movie1.time_stamps(end)], [handles.recording.abf.Wheel(b),handles.recording.abf.Wheel(b)],'Color','b');
line([handles.recording.movie1.time_stamps(1),handles.recording.movie1.time_stamps(end)], [otherbound,otherbound],'Color','g');

axes(handles.axes5);
line([handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value)),handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))], [-10,10],'Color','m');


end


% --- Executes on slider movement.
function slidermin_Callback(hObject, eventdata, handles)
% hObject    handle to slidermin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function slidermin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidermin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function slidermax_Callback(hObject, eventdata, handles)
% hObject    handle to slidermax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function slidermax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidermax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on selection change in cmap_select2.
function cmap_select2_Callback(hObject, eventdata, handles)
% hObject    handle to cmap_select2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap_select2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap_select2
populate_graphs(hObject);

end

% --- Executes during object creation, after setting all properties.
function cmap_select2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmap_select2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in sucrose_roi.
function sucrose_roi_Callback(hObject, eventdata, handles)
% hObject    handle to sucrose_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sucrose_roi

plot_sucrose_ROI(hObject);

end


% --- Executes on slider movement.
function sucrose_roi_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to sucrose_roi_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
plot_sucrose_ROI(hObject);

end

% --- Executes during object creation, after setting all properties.
function sucrose_roi_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sucrose_roi_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes on button press in shrinkaxis.
function shrinkaxis_Callback(hObject, eventdata, handles)
% hObject    handle to shrinkaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shrinkaxis
axes(handles.axes3);  hold on;

if(handles.shrinkaxis.Value == 0)
    set(gca,'xlim',[handles.time_st, handles.time_en]);
else
    set(gca,'xlim',[handles.recording.movie1.time_stamps(1), handles.recording.movie1.time_stamps(end)]);
end

axes(handles.axes5);  hold on;

if(handles.shrinkaxis.Value == 0)
    set(gca,'xlim',[handles.time_st, handles.time_en]);
else
    set(gca,'xlim',[handles.recording.movie1.time_stamps(1), handles.recording.movie1.time_stamps(end)]);
end

end


% --- Executes on button press in newroi.
function newroi_Callback(hObject, eventdata, handles)
% hObject    handle to newroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); hold on;
hFH = imfreehand(gca);
binaryImage = hFH.createMask();
handles.total_ROI = handles.total_ROI + 1;
set(handles.total_roi_display, 'String', num2str(handles.total_ROI));

handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) = binaryImage;

guidata(hObject, handles);
end

% --- Executes on button press in addtoroi.
function addtoroi_Callback(hObject, eventdata, handles)
% hObject    handle to addtoroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); hold on;

hFH = imfreehand(gca);
binaryImage = hFH.createMask();
set(handles.total_roi_display, 'String', num2str(handles.total_ROI));

tmpmat = handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) + binaryImage;
tmpmat(tmpmat > 1) = 1;
handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) = tmpmat;

guidata(hObject, handles);
end

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% calculate df over f etc in different ROIs
size_tsStack = size(handles.recording.tseries(handles.tseries_for_plot).tsStack);
tsStack = handles.recording.tseries(handles.tseries_for_plot).tsStack;

for k =1:1:handles.total_ROI
    for j = 1:1:handles.recording.tseries(handles.tseries_for_plot).Total2pFrames
        if(length(size_tsStack) == 5)
            tmp = tsStack(:,:,1,j,:);
            tmp = squeeze(tmp);
            tmp_masked_data = tmp(handles.ROI_masks(:,:,:,k) == 1);
            [a b]  = find(handles.ROI_masks(:,:,:,k) == 1);
            handles.custom_ROI_mean(k,j) = sum(tmp_masked_data)./length(a);
            
        else
            tmp = tsStack(:,:,j);
            tmp = squeeze(tmp);
            tmp_masked_data = tmp(handles.ROI_masks(:,:,1,k) == 1);
            [a b]  = find(handles.ROI_masks(:,:,1,k) == 1);
            handles.custom_ROI_mean(k,j) = sum(tmp_masked_data)./length(a);
        end
    end
end

for  k =1:1:handles.total_ROI
    handles.custom_ROI_df_over_f_usingmean(k,:) = (handles.custom_ROI_mean(k,:)-mean(handles.custom_ROI_mean(k,:)))./mean(handles.custom_ROI_mean(k,:));
    low_val = prctile( handles.custom_ROI_mean(k,:),5);
    handles.custom_ROI_df_over_f(k,:) = (handles.custom_ROI_mean(k,:)-low_val)./low_val;
end

%plot new ROI from 2p
axes(handles.axes3); cla('reset'); hold on;
imagesc(handles.recording.tseries(handles.tseries_for_plot).Time_s,1:handles.total_ROI,handles.custom_ROI_df_over_f);
ylabel('ROIs','fontsize', 8);
set(gca,'ylim',[0.5 handles.total_ROI+.5]);
set(gca,'xlim',[handles.time_st, handles.time_en]);
set(gca,'TickDir','out')
set(gca,'TickLength',[0.004 0.004]);
set(gca,'FontSize',8)
box on;
line([handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value)),handles.recording.movie1.time_stamps(ceil(handles.t_slider.Value))], [-10,10],'Color','m');


end

% --- Executes on button press in deleteroi.
function deleteroi_Callback(hObject, eventdata, handles)
% hObject    handle to deleteroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.total_ROI > 0)
    handles.ROI_masks(:,:,:,handles.total_ROI) = 0;
    handles.total_ROI = handles.total_ROI - 1;
    set(handles.total_roi_display, 'String', num2str(handles.total_ROI));
end
guidata(hObject, handles);
end


% --- Executes on button press in plotselectedrois.
function plotselectedrois_Callback(hObject, eventdata, handles)
% hObject    handle to plotselectedrois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotselectedrois
end


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global play_movie;
play_movie = 1;

% if(handles.stop_play)
%     handles.stop_play = 0;
% else
%     handles.stop_play = 1;
% end

while(play_movie)
    set(handles.t_slider, 'Value', min([handles.t_slider.Value+str2num(handles.fps.String),length(handles.recording.movie1.time_stamps)]) );
    set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
    populate_graphs(hObject);
    pause(.1);
end

end



function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double
end

% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.stop, 'Value', ~ get(set(handles.stop, 'Value')))
global play_movie;
global MOVIE;

close(MOVIE);
play_movie = 0;


% if(handles.stop_play)
%     handles.stop_play = 0;
% else
%     handles.stop_play = 1;
% end
%guidata(hObject, handles);

end


% --- Executes on button press in mean_z.
function mean_z_Callback(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mean_z
end


% --- Executes on button press in record.
function record_Callback(hObject, eventdata, handles)
% hObject    handle to record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global play_movie;
global MOVIE;
handles.recording.abfname
MOVIE = VideoWriter([handles.recording.abfname(1:end-4) '_movie.avi']);
MOVIE.FrameRate = 25./str2num(handles.fps.String);
open(MOVIE);
play_movie = 1;

% if(handles.stop_play)
%     handles.stop_play = 0;
% else
%     handles.stop_play = 1;
% end

while(play_movie)
    set(handles.t_slider, 'Value', min([handles.t_slider.Value+str2num(handles.fps.String),length(handles.recording.movie1.time_stamps)]) );
    set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
    populate_graphs(hObject);
    %pause(.1);
    drawnow;
    frame = getframe(gcf);
    writeVideo(MOVIE,frame);
    
end


end
