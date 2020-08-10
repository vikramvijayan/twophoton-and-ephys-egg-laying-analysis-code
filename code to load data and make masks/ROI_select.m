function varargout = ROI_select(varargin)
% ROI_SELECT MATLAB code for ROI_select.fig
%      ROI_SELECT, by itself, creates a new ROI_SELECT or raises the existing
%      singleton*.
%
%      H = ROI_SELECT returns the handle to a new ROI_SELECT or the handle to
%      the existing singleton*.
%
%      ROI_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_SELECT.M with the given input arguments.
%
%      ROI_SELECT('Property','Value',...) creates a new ROI_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROI_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROI_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_select

% Last Modified by GUIDE v2.5 08-May-2019 13:34:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ROI_select_OpeningFcn, ...
    'gui_OutputFcn',  @ROI_select_OutputFcn, ...
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
% End initialization code - DO NOT EDIT
warning('off');

end

% --- Executes just before ROI_select is made visible.
function ROI_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROI_select (see VARARGIN)


InputStack = varargin{1};
FileName = varargin{2};
TakeMax = varargin{3};

%if TakeMax, then we will be taking the ROI from the first stack and using
%it for all stacks

handles.TakeMax = TakeMax;
handles.tsStack = InputStack;
handles.filename = FileName;

size_mat = size(InputStack);

handles.size_mat = size_mat;



handles.total_ROI = 0;

% if there is more than 1 z slice
if(length(size_mat) == 5)
    handles.ROI_masks = zeros(size_mat(1), size_mat(2), size_mat(5),20);

    
    set(handles.z_slider,'Max',handles.size_mat(5),'Min',1);
    set(handles.z_slider,'SliderStep',[1/(size_mat(5)-1) 1/(size_mat(5)-1)]);
    set(handles.z_slider, 'Value', 1);
    
    set(handles.t_slider,'Max',handles.size_mat(4),'Min',1);
    set(handles.t_slider,'SliderStep',[1/(size_mat(4)-1) 10/(size_mat(4)-1)]);
    set(handles.t_slider, 'Value', 1);
    
else
    
    handles.ROI_masks = zeros(size_mat(1), size_mat(2), 1,20);
    
    
    set(handles.z_slider,'Max',1,'Min',1);
    set(handles.z_slider,'SliderStep',[0 0]);
    set(handles.z_slider, 'Value', 1);
    
    set(handles.t_slider,'Max',handles.size_mat(3),'Min',1);
    set(handles.t_slider,'SliderStep',[1/(size_mat(3)-1) 10/(size_mat(3)-1)]);
    set(handles.t_slider, 'Value', 1);
    
    
end

set(handles.current_t_display, 'String', num2str(ceil(handles.t_slider.Value)));
set(handles.current_z_display, 'String', num2str(ceil(handles.z_slider.Value)));
set(handles.total_roi_display, 'String', num2str(handles.total_ROI));

set(handles.maxcolor,'Max',2^10,'Min',0);
set(handles.maxcolor,'SliderStep',[1/(2^10-1) 10/(2^10-1)]);
set(handles.maxcolor, 'Value', 2^10);

set(handles.mincolor,'Max',2^10,'Min',0);
set(handles.mincolor,'SliderStep',[1/(2^10-1) 10/(2^10-1)]);
set(handles.mincolor, 'Value', 0);

if(max(max(max(max(max((handles.tsStack(:,:,:,:,:))))))) > 2^12)
    set(handles.maxcolor,'Max',2^16,'Min',0);
    set(handles.maxcolor,'SliderStep',[1/(2^16-1) 10/(2^16-1)]);
    set(handles.maxcolor, 'Value', 2^16);
    
    set(handles.mincolor,'Max',2^16,'Min',0);
    set(handles.mincolor,'SliderStep',[1/(2^16-1) 10/(2^16-1)]);
    set(handles.mincolor, 'Value', 0);
end
set(handles.cmap_select, 'String', {'parula','jet','hsv','hot','cool','gray'});
handles.cmap_select.Value = 1;
handles.cmap_options  =  {'parula','jet','hsv','hot','cool','gray'};

% Choose default command line output for ROI_select
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

plot_ROI(handles);

% UIWAIT makes ROI_select wait for user response (see UIRESUME)
end

% --- Outputs from this function are returned to the command line.
function varargout = ROI_select_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
warning('on');


varargout{1} = handles.ROI_masks;
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
plot_ROI(handles);

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
plot_ROI(handles);

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


% --- Executes on button press in new_roi.
function new_roi_Callback(hObject, eventdata, handles)
% hObject    handle to new_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hFH = imfreehand();
binaryImage = hFH.createMask();
handles.total_ROI = handles.total_ROI + 1;
set(handles.total_roi_display, 'String', num2str(handles.total_ROI));

handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) = binaryImage;

guidata(hObject, handles);
%plot_ROI(handles);

end
% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.ROI_masks;
ROI_masks_from_GUI = handles.ROI_masks(:,:,:,1:handles.total_ROI);

if(handles.TakeMax)
    for i = 1:1:handles.total_ROI
        for j = 2:1:handles.size_mat(5)
            ROI_masks_from_GUI(:,:,j,i) = ROI_masks_from_GUI(:,:,1,i);
        end
    end
end

save([handles.filename '.mat'], 'ROI_masks_from_GUI');
close(handles.figure1)
end


% --- Executes on button press in delete_roi.
function delete_roi_Callback(hObject, eventdata, handles)
% hObject    handle to delete_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.total_ROI > 0)
    handles.ROI_masks(:,:,:,handles.total_ROI) = 0;
    handles.total_ROI = handles.total_ROI - 1;
    set(handles.total_roi_display, 'String', num2str(handles.total_ROI));
end
guidata(hObject, handles);
%plot_ROI(handles);

end

function plot_ROI(handles)

map = colormap(char(handles.cmap_options(handles.cmap_select.Value)));
set(gca, 'CLim', [1, 2^16-1]);
colormap(map);

% if there are z
if(length(handles.size_mat) == 5)
    
    if(~handles.max_t.Value && ~handles.max_z.Value)
        curr_image = handles.tsStack(:,:,end,round(handles.t_slider.Value), round(handles.z_slider.Value));
    end
    
    if(handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,end,:,round(handles.z_slider.Value)),[],4);
    end
    
    if(handles.max_z.Value)
        curr_image = max(handles.tsStack(:,:,end,round(handles.t_slider.Value),:),[],5);
    end
    
    if(handles.mean_t.Value)
        curr_image = mean(handles.tsStack(:,:,end,:,round(handles.z_slider.Value)),4);
    end
    
    if(handles.mean_z.Value)
        curr_image = mean(handles.tsStack(:,:,end,round(handles.t_slider.Value),:),5);
    end
    
    if(handles.mean_z.Value & handles.mean_t.Value)
        curr_image = mean(mean(handles.tsStack(:,:,end,:,:),4),5);
    end
    
    if(handles.max_z.Value & handles.max_t.Value)
        curr_image = max(max(handles.tsStack(:,:,end,:,:),[],5),[],4);
    end
    
    fl = 0;
    
    
    
    if(handles.show_roi.Value & handles.total_ROI > 0 & handles.max_t.Value  & handles.max_z.Value)
        for i = 1:1:handles.total_ROI
            curr_image(max(handles.ROI_masks(:,:,:,i),[],3) == 1) = ceil((2^16-1)/i);
        end
        imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
        fl = 1;
    end
    
    if(handles.show_roi.Value & handles.total_ROI > 0  & ~handles.max_t.Value & handles.max_z.Value)
        for i = 1:1:handles.total_ROI
            curr_image(max(handles.ROI_masks(:,:,:,i),[],3) == 1) = ceil((2^16-1)/i);
        end
        imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
        fl = 1;
    end
    
    if(handles.show_roi.Value & handles.total_ROI > 0  & ~handles.max_z.Value)
        for i = 1:1:handles.total_ROI
            curr_image(handles.ROI_masks(:,:,handles.z_slider.Value,i) == 1) = ceil((2^16-1)/i);
        end
        imshow(curr_image(:,:),[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
        fl = 1;
    end
    
    
    
    if(~fl)
        imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
    end
    
else
      if(~handles.max_t.Value && ~handles.max_z.Value)
        curr_image = handles.tsStack(:,:,round(handles.t_slider.Value));
    end
    
    if(handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,:),[],3);
    end
    
    if(handles.max_z.Value)
        curr_image = handles.tsStack(:,:,round(handles.t_slider.Value));
    end
    
    if(handles.max_z.Value & handles.max_t.Value)
        curr_image = max(handles.tsStack(:,:,:),[],3);
    end
    
       if(handles.mean_t.Value)
        curr_image = mean(handles.tsStack(:,:,:),3);
    end
    
    
    fl = 0;
    
    
    
    if(handles.show_roi.Value & handles.total_ROI > 0)
        for i = 1:1:handles.total_ROI
            curr_image(handles.ROI_masks(:,:,1,1) == 1) = ceil((2^16-1)/i);
        end
        imshow(curr_image(:,:),[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
        fl = 1;
    end
    

    if(~fl)
        imshow(curr_image,[round(handles.mincolor.Value), round(handles.maxcolor.Value)]);
    end
    

end

colormap(map);

end


% --- Executes on button press in show_roi.
function show_roi_Callback(hObject, eventdata, handles)
% hObject    handle to show_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_ROI(handles);
% Hint: get(hObject,'Value') returns toggle state of show_roi
end


% --- Executes on button press in max_t.
function max_t_Callback(hObject, eventdata, handles)
% hObject    handle to max_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_ROI(handles);

% Hint: get(hObject,'Value') returns toggle state of max_t
end

% --- Executes on button press in max_z.
function max_z_Callback(hObject, eventdata, handles)
% hObject    handle to max_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_ROI(handles);

% Hint: get(hObject,'Value') returns toggle state of max_z
end


% --- Executes on selection change in cmap_select.
function cmap_select_Callback(hObject, eventdata, handles)
% hObject    handle to cmap_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap_select
plot_ROI(handles)
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


% --- Executes on button press in addtoroi.
function addtoroi_Callback(hObject, eventdata, handles)
% hObject    handle to addtoroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hFH = imfreehand();
binaryImage = hFH.createMask();
set(handles.total_roi_display, 'String', num2str(handles.total_ROI));

tmpmat = handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) + binaryImage;
tmpmat(tmpmat > 1) = 1;
handles.ROI_masks(:,:,round(handles.z_slider.Value),handles.total_ROI) = tmpmat;

guidata(hObject, handles);
%plot_ROI(handles);
end


% --- Executes on slider movement.
function mincolor_Callback(hObject, eventdata, handles)
% hObject    handle to mincolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
guidata(hObject, handles);
plot_ROI(handles);
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

guidata(hObject, handles);
plot_ROI(handles);
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


% --- Executes on button press in mean_t.
function mean_t_Callback(hObject, eventdata, handles)
% hObject    handle to mean_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mean_t
end

% --- Executes on button press in mean_z.
function mean_z_Callback(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mean_z
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = mean(handles.tsStack(:,:,end,:,:),5);
tmp2 = mean(tmp,4);
tmp2 = squeeze(tmp2);

automask(tmp2,handles.filename,handles.size_mat)
end
