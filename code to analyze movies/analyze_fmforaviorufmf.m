function varargout = analyze_fmforaviorufmf(varargin)
% ANALYZE_FMFORAVIORUFMF MATLAB code for analyze_fmforaviorufmf.fig
%      ANALYZE_FMFORAVIORUFMF, by itself, creates a new ANALYZE_FMFORAVIORUFMF or raises the existing
%      singleton*.
%
%      H = ANALYZE_FMFORAVIORUFMF returns the handle to a new ANALYZE_FMFORAVIORUFMF or the handle to
%      the existing singleton*.
%
%      ANALYZE_FMFORAVIORUFMF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZE_FMFORAVIORUFMF.M with the given input arguments.
%
%      ANALYZE_FMFORAVIORUFMF('Property','Value',...) creates a new ANALYZE_FMFORAVIORUFMF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyze_fmforaviorufmf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyze_fmforaviorufmf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyze_fmforaviorufmf

% Last Modified by GUIDE v2.5 06-Feb-2018 11:58:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @analyze_fmforaviorufmf_OpeningFcn, ...
    'gui_OutputFcn',  @analyze_fmforaviorufmf_OutputFcn, ...
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

% --- Executes just before analyze_fmforaviorufmf is made visible.
function analyze_fmforaviorufmf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyze_fmforaviorufmf (see VARARGIN)


warning('off');

handles.centerplainval = -1;
handles.filename = varargin{1};

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
    handles.max_n_frames = handles.vobj.FrameRate*handles.vobj.Duration;
    handles.movie_length = handles.vobj.FrameRate*handles.vobj.Duration;
    handles.f_height = handles.vobj.Height;
    handles.f_width = handles.vobj.Width;
    handles.isfmf = 0;
    
end

handles.binaryImage = zeros(handles.f_height,handles.f_width);
handles.egg_array = {};
set(handles.listbox1, 'String', handles.egg_array);

%handles.first_stamp = stamp;

set(handles.t_slider,'Max',handles.movie_length,'Min',1);
set(handles.t_slider,'SliderStep',[1/(handles.movie_length-1) 18/(handles.movie_length-1)]);
set(handles.t_slider, 'Value', 1);

set(handles.cmap_select, 'String', {'gray','parula','jet','hsv','hot','cool'});
handles.cmap_select.Value = 1;
handles.cmap_options  =  {'gray','parula','jet','hsv','hot','cool','gray'};

set(handles.slidermax,'Max',2^8,'Min',0);
set(handles.slidermax,'SliderStep',[1/(2^8-1) 10/(2^8-1)]);
set(handles.slidermax, 'Value', 2^8);

set(handles.slidermin,'Max',2^8,'Min',0);
set(handles.slidermin,'SliderStep',[1/(2^8-1) 10/(2^8-1)]);
set(handles.slidermin, 'Value', 0);


set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));


handles.first_run = 1;

% was having some problems with scroll wheel so disabled for the moment
% set(gcf, 'WindowScrollWheelFcn', {@wheel,handles});
% handles.sliderListener = addlistener(handles.t_slider,'ContinuousValueChange', @(hObject, event) t_slider_Callback(hObject, eventdata, handles));

% Choose default command line output for view_tsstack_egg
handles.output = hObject;

% save the changes to handles
guidata(hObject, handles);

% UIWAIT makes analyze_fmforaviorufmf wait for user response (see UIRESUME)
end

% --- Outputs from this function are returned to the command line.
function varargout = analyze_fmforaviorufmf_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
warning('on');

%delete(handles.figure1);
end



% --- Executes on slider movement.
function t_slider_Callback(hObject, eventdata, handles)
% hObject    handle to t_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
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

if(handles.isfmf == 1)
    plot_movie_fmf(hObject);
end


if(handles.isfmf == 2)
    plot_movie_ufmf(hObject);
end


if(handles.isfmf == 0)
    plot_movie_avi(hObject);
end


end



function plot_movie_fmf(hObject)
handles = guidata(hObject);


fseek(handles.fp,handles.header_size+(ceil(handles.t_slider.Value)-1)*handles.bytes_per_chunk,'bof');
[data, stamp] = fmf_read_frame( handles.fp, handles.f_height, handles.f_width, handles.bytes_per_chunk, handles.data_format );

% need to change this to the stamp correcte for abf start
% set(handles.moviestamp, 'String', num2str(stamp-handles.first_stamp));

vidFrame = double(data) + handles.binaryImage.*2^8;

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
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
else
    axes(handles.axes2); hold on;
    set(handles.movie_show ,'cdata',flipud(vidFrame));
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end

guidata(hObject, handles);
end


function plot_movie_ufmf(hObject)
handles = guidata(hObject);

[im,header,timestamp,bb,mu] = ufmf_read_frame(handles.header,ceil(handles.t_slider.Value));


% need to change this to the stamp correcte for abf start
% set(handles.moviestamp, 'String', num2str(stamp-handles.first_stamp));

vidFrame = double(im) + handles.binaryImage.*2^8;

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
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
else
    axes(handles.axes2); hold on;
    set(handles.movie_show ,'cdata',flipud(vidFrame));
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end

guidata(hObject, handles);
end




function plot_movie_avi(hObject)
handles = guidata(hObject);


handles.vobj.CurrentTime = ceil(handles.t_slider.Value)./handles.vobj.FrameRate - 1/handles.vobj.FrameRate;
vidFrame3 = readFrame(handles.vobj);

vidFrame = double(vidFrame3(:,:,1)) + handles.binaryImage.*2^8;



if(handles.first_run)
    axes(handles.axes2); hold on;
    handles.movie_show = imagesc(flipud(vidFrame));
    axis tight;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    handles.first_run = 0;
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
else
    axes(handles.axes2); hold on;
    set(handles.movie_show ,'cdata',flipud(vidFrame));
    map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
    colormap(map);
    set(gca, 'CLim', [handles.slidermin.Value, handles.slidermax.Value]);
end


guidata(hObject, handles);
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


% --- Executes on button press in sucrose.
function sucrose_Callback(hObject, eventdata, handles)
% hObject    handle to sucrose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
hFH = imfreehand();
binaryImage = hFH.createMask();


binaryImage = flipud(binaryImage);


handles.binaryImage = binaryImage;

guidata(hObject, handles);
%populate_graphs(handles);


end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stamp = zeros(1,uint16(handles.max_n_frames));
roi_sucrose = zeros(1,uint16(handles.max_n_frames));

if(handles.isfmf == 1)
    fseek(handles.fp,handles.header_size,'bof');
    
    for cnt = 1:1:handles.max_n_frames
        [dataout, stampout] = fmf_read_frame( handles.fp, handles.f_height, handles.f_width, handles.bytes_per_chunk, handles.data_format );
        stampfmf(cnt) = stampout;
        roi_sucrose(cnt) = sum(sum(handles.binaryImage.*double(dataout)));
        roi_body2(cnt) = sum(sum(handles.binaryImage2.*double(dataout)));
        roi_body1(cnt) = sum(sum(handles.binaryImage3.*double(dataout)));
        
        m(:) = dataout(1,1:4,1);
        mm(:) = convert_d(double(m(:)));
        cyclecount(cnt) = bin2dec(mm(8:20));
        secondcount(cnt) = bin2dec(mm(1:7));
        cycleoff(cnt) = bin2dec(mm(21:32));
    end
    stamps = secondcount + (cyclecount+cycleoff./3072)./8000;
    
    add = 0;
    addnum = 0;
    
    for i = 1:1:(length(stamps)-1)
        if(stamps(i+1)-stamps(i) < 0)
            addnum = addnum+128;
        end
        add = [add addnum];
    end
    
    stamp = stamps+add;
end

if(handles.isfmf == 2)
    
    for cnt = 1:1:handles.max_n_frames
        [im,header,stampout,bb,mu] = ufmf_read_frame(handles.header,cnt);
        stampfmf(cnt) = stampout;
        roi_sucrose(cnt) = sum(sum(handles.binaryImage.*double(im)));
        roi_body1(cnt) = sum(sum(handles.binaryImage2.*double(im)));
        roi_body2(cnt) = sum(sum(handles.binaryImage3.*double(im)));
        
        m(:) = im(1,1:4,1);
        mm(:) = convert_d(double(m(:)));
        cyclecount(cnt) = bin2dec(mm(8:20));
        secondcount(cnt) = bin2dec(mm(1:7));
        cycleoff(cnt) = bin2dec(mm(21:32));
    end
    stamps = secondcount + (cyclecount+cycleoff./3072)./8000;
    
    add = 0;
    addnum = 0;
    
    for i = 1:1:(length(stamps)-1)
        if(stamps(i+1)-stamps(i) < 0)
            addnum = addnum+128;
        end
        add = [add addnum];
    end
    
    stamp = stamps+add;
end

if(handles.isfmf == 0)
    stampfmf = [];
    handles.vobj.CurrentTime = 0;
    len = 0;
    while(hasFrame(handles.vobj))
        len = len+1;
        q(:,:,:) = readFrame(handles.vobj);
        m(:) = q(1,1:4,1);
        mm(:) = convert_d(double(m(:)));
        cyclecount(len) = bin2dec(mm(8:20));
        secondcount(len) = bin2dec(mm(1:7));
        cycleoff(len) = bin2dec(mm(21:32));
        roi_sucrose(len) = sum(sum(handles.binaryImage.*double(q(:,:,1))));
        roi_body1(len) = sum(sum(handles.binaryImage2.*double(q(:,:,1))));
        roi_body2(len) = sum(sum(handles.binaryImage3.*double(q(:,:,1))));
        
    end
    
    stamps = secondcount + (cyclecount+cycleoff./3072)./8000;
    
    add = 0;
    addnum = 0;
    
    for i = 1:1:(length(stamps)-1)
        if(stamps(i+1)-stamps(i) < 0)
            addnum = addnum+128;
        end
        add = [add addnum];
    end
    
    stamp = stamps+add;
end


egg_array_num = [];
for i = 1:1:length(handles.egg_array)
    egg_array_num(i) = (cell2mat((handles.egg_array(i))));
end

egg_array_num = egg_array_num';
stamp = stamp';
stampfmf = stampfmf';
roi_sucrose = roi_sucrose';
roi_body1 = roi_body1';
roi_body2 = roi_body2';

bound1 = handles.bound1val;
bound2 = handles.bound2val;
marksucrose = handles.marksucroseval;

centerplain = handles.centerplainval;

modifiedStr = strrep([handles.filename '_processed'], '.', '_');
save(modifiedStr,'stamp','roi_sucrose','egg_array_num','stampfmf','bound1','bound2','marksucrose','roi_body1','roi_body2','centerplain');

set(handles.saving, 'String', 'Done Saving');

end

% --- Executes on button press in add_egg.
function add_egg_Callback(hObject, eventdata, handles)
% hObject    handle to add_egg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.egg_array (end+1)= num2cell(ceil(handles.t_slider.Value));
set(handles.listbox1, 'String', handles.egg_array);
guidata(hObject, handles);

end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

contents = cellstr(get(hObject,'String'));
current_select = contents{get(hObject,'Value')};
set(handles.t_slider, 'Value', (str2num((current_select))));
set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
guidata(hObject, handles);
populate_graphs(hObject);


end

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in rem_egg.
function rem_egg_Callback(hObject, eventdata, handles)
% hObject    handle to rem_egg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.egg_array = handles.egg_array(1:end-1);
set(handles.listbox1, 'String', handles.egg_array);
guidata(hObject, handles);

end


% --- Executes on slider movement.
function slidermax_Callback(hObject, eventdata, handles)
% hObject    handle to slidermax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
guidata(hObject, handles);

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
% --- Executes on slider movement.
function slidermin_Callback(hObject, eventdata, handles)
% hObject    handle to slidermin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
guidata(hObject, handles);

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


% --- Executes on selection change in cmap_select.
function cmap_select_Callback(hObject, eventdata, handles)
% hObject    handle to cmap_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap_select
guidata(hObject, handles);

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


% --- Executes on button press in bound1.
function bound1_Callback(hObject, eventdata, handles)
% hObject    handle to bound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bound1val = ceil(handles.t_slider.Value);
set(handles.bound1text,'String',handles.bound1val);
guidata(hObject, handles);

end

% --- Executes on button press in bound2.
function bound2_Callback(hObject, eventdata, handles)
% hObject    handle to bound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bound2val = ceil(handles.t_slider.Value);
set(handles.bound2text,'String',ceil(handles.bound2val));
guidata(hObject, handles);

end

% --- Executes on button press in marksucrose.
function marksucrose_Callback(hObject, eventdata, handles)
% hObject    handle to marksucrose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.marksucroseval = ceil(handles.t_slider.Value);
set(handles.marksucrosetext,'String',ceil(handles.marksucroseval));
guidata(hObject, handles);

end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
hFH = imfreehand();
binaryImage = hFH.createMask();


binaryImage = flipud(binaryImage);


handles.binaryImage2 = binaryImage;

guidata(hObject, handles);
end

% --- Executes on button press in bodyROI2. (actually body 1)
function bodyROI2_Callback(hObject, eventdata, handles)
% hObject    handle to bodyROI2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
hFH = imfreehand();
binaryImage = hFH.createMask();


binaryImage = flipud(binaryImage);


handles.binaryImage3 = binaryImage;

guidata(hObject, handles);
end


% --- Executes on button press in centerplain.
function centerplain_Callback(hObject, eventdata, handles)
% hObject    handle to centerplain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.centerplainval = ceil(handles.t_slider.Value);
set(handles.centerplaintext,'String',ceil(handles.centerplainval));
guidata(hObject, handles);

end
