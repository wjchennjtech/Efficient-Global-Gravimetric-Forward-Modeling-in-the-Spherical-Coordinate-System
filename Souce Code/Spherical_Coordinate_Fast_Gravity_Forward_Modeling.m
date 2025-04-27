function varargout = Spherical_Coordinate_Fast_Gravity_Forward_Modeling(varargin)
% SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING MATLAB code for Spherical_Coordinate_Fast_Gravity_Forward_Modeling.fig
%      SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING, by itself, creates a new SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING or raises the existing
%      singleton*.
%
%      H = SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING returns the handle to a new SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING or the handle to
%      the existing singleton*.
%
%      SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING.M with the given input arguments.
%
%      SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING('Property','Value',...) creates a new SPHERICAL_COORDINATE_FAST_GRAVITY_FORWARD_MODELING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spherical_Coordinate_Fast_Gravity_Forward_Modeling

% Last Modified by GUIDE v2.5 18-Oct-2023 08:59:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OpeningFcn, ...
                   'gui_OutputFcn',  @Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OutputFcn, ...
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


% --- Executes just before Spherical_Coordinate_Fast_Gravity_Forward_Modeling is made visible.
function Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spherical_Coordinate_Fast_Gravity_Forward_Modeling (see VARARGIN)

% Choose default command line output for Spherical_Coordinate_Fast_Gravity_Forward_Modeling
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Spherical_Coordinate_Fast_Gravity_Forward_Modeling wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.edit9, 'enable', 'off');
set(handles.edit9, 'string', ' ');
set(handles.pushbutton14,'enable','of');

% --- Outputs from this function are returned to the command line.
function varargout = Spherical_Coordinate_Fast_Gravity_Forward_Modeling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

set(handles.edit7, 'enable', 'on');
set(handles.edit8, 'enable', 'on');
set(handles.edit9, 'enable', 'off');
set(handles.edit9, 'string', ' ');

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

set(handles.edit7, 'enable', 'on');
set(handles.edit8, 'enable', 'on');
set(handles.edit9, 'enable', 'on');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data;global nlon;global nlat;global nlayers;


nmax=str2double(get(handles.edit7,'string'));
caculat_height=str2double(get(handles.edit8,'string'));
Lat=0;
Lon=0;
if handles.radiobutton6.Value == 1
lon_min=str2double(get(handles.edit1,'string'));
lon_max=str2double(get(handles.edit2,'string'));
lat_min=str2double(get(handles.edit4,'string'));
lat_max=str2double(get(handles.edit5,'string'));

dlon=abs(str2double(get(handles.edit3,'string')));
dlat=abs(str2double(get(handles.edit6,'string')));
if isnan(dlat)
dlat=abs(eval(get(handles.edit6,'string')));
end
if isnan(dlon)
dlon=abs(eval(get(handles.edit3,'string')));
end
nlon=360/dlon;
nlat=180/dlat;    
Lat=linspace(lat_max,lat_min,nlat);
Lon=linspace(lon_min,lon_max,nlon);

elseif handles.radiobutton7.Value == 1
data_table=get(handles.uitable7,'data');
Lat=data_table(:,2);
Lon=data_table(:,1);
end

%%
%estimate the input data
lonmax=max(data(:,1));lonmin=min(data(:,1));
latmax=max(data(:,2));latmin=min(data(:,2));
ifzero=0;
for i=1:size(data,1)
Dlon1=rem(abs(data(1,1)-data(2,1)),1);
Dlon2=rem(abs(data(1,1)-data(i+2,1)),1);
if Dlon1==0 && Dlon2==0
Dlon=abs(data(1,1)-data(2,1));
if any(rem(lonmin,1))
Nlon=360/Dlon;
else
Nlon=360/Dlon+1;
ifzero=1;
end
break;
end

if Dlon1>Dlon2
numb_lon=i+2;
Dlon=1/(numb_lon-1);

if any(rem(lonmin,1))
Nlon=360/Dlon;
else
Nlon=360/Dlon+1;
ifzero=1;
end
break;
end
end

Nlat=size(data,1)/Nlon/nlayers;
if ifzero==0
Dlat=180/Nlat;
else
Dlat=180/(Nlat-1);    
end
%%
f=[];
if handles.radiobutton6.Value == 1
fresult=zeros(nlat*nlon,3);
elseif  handles.radiobutton7.Value == 1
fresult=zeros(size(Lat,1),3);
end

gravity_potential=0;
gravity_attraction=0;
gravity_gradient=0;
if handles.radiobutton1.Value == 1
%straight
nrow=size(data,1)/nlayers;
for i=1:nlayers
%to computa multi-layers 
DATA=data((i-1)*nrow+1:i*nrow,:);
[result_cnm_snm]=analysis_straight(nmax,Dlat,Dlon,Nlat,Nlon,DATA);
if handles.radiobutton6.Value == 1
[Gravity_potential,Gravity_attraction,Gravity_gradient]=synthesis_straight(result_cnm_snm,nmax,lat_max,lat_min,lon_max,lon_min,dlat,dlon,caculat_height);
elseif handles.radiobutton7.Value == 1
[Gravity_potential,Gravity_attraction,Gravity_gradient]=synthesis_straight_point(result_cnm_snm,nmax,Lat,Lon,caculat_height);    
end
gravity_potential=gravity_potential+Gravity_potential;
gravity_attraction=gravity_attraction+Gravity_attraction;
gravity_gradient=gravity_gradient+Gravity_gradient;
end
if handles.radiobutton3.Value == 1

f=gravity_potential;

else if handles.radiobutton4.Value == 1

f=gravity_attraction*1e+5;

else if handles.radiobutton5.Value == 1

f=gravity_gradient*1e+9;

end
end
end

else if  handles.radiobutton2.Value == 1
%slice

num_latyers=str2double(get(handles.edit9,'string'));
[Cnm,Snm]=analysis_slice(nmax,Dlat,Dlon,Nlat,Nlon,data,num_latyers);

if handles.radiobutton6.Value == 1

[gravity_potential,gravity_attraction,gravity_gradient]=synthesis_slice(Cnm,Snm,nmax,lat_max,lat_min,lon_max,lon_min,dlat,dlon,caculat_height,num_latyers);

elseif handles.radiobutton7.Value == 1
   
[gravity_potential,gravity_attraction,gravity_gradient]=synthesis_slice_point(Cnm,Snm,nmax,Lat,Lon,caculat_height,num_latyers);
end
if handles.radiobutton3.Value == 1

f=gravity_potential;

else if handles.radiobutton4.Value == 1

f=gravity_attraction*1e+5;

else if handles.radiobutton5.Value == 1

f=gravity_gradient*1e+9;

end
end
end

end
end
if handles.radiobutton6.Value == 1

    for i=1:nlat
        for j=1:nlon

            fresult((i-1)*nlon+j,1)=Lon(j);
            fresult((i-1)*nlon+j,2)=Lat(i);
            fresult((i-1)*nlon+j,3)=f(i,j);

        end
    end
elseif handles.radiobutton7.Value == 1

fresult(:,1)=Lon;
fresult(:,2)=Lat;
fresult(:,3)=f;

end
fmin=min(fresult(:,3));
fmax=max(fresult(:,3));
fmean=mean(fresult(:,3));
fstd=std(fresult(:,3));
stastic=[fmin,fmax,fmean,fstd];
set(handles.uitable2,'data',stastic);

set(handles.uitable1,'data',fresult);    



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data1=Model_data;
set(data1,'visible','on');

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nlon;global nlat;
table1_data = get(handles.uitable1,'data');
if  handles.radiobutton6.Value == 1
latmax=ceil(max(table1_data(:,2)));
latmin=floor(min(table1_data(:,2)));
lonmax=ceil(max(table1_data(:,1)));
lonmin=floor(min(table1_data(:,1)));

datap=reshape(table1_data(:,3),nlon,nlat);
data_p=datap';
data_p=flipud(data_p);
imagesc([lonmin lonmax],[latmin latmax],data_p);

xlim([lonmin lonmax]);
ylim([latmin latmax]);

xticks(lonmin:90:lonmax);
yticks(latmin:45:latmax);
 set(gca, 'YDir', 'normal');
elseif  handles.radiobutton7.Value == 1
datap=table1_data;
lonmax=ceil(max(datap(:,1)));
lonmin=floor(min(datap(:,1)));
latmax=ceil(max(datap(:,2)));
latmin=floor(min(datap(:,2)));

scatter(datap(:, 1), datap(:, 2), 1, datap(:, 3));

xlim([lonmin lonmax]);
ylim([latmin latmax]);

xticks(lonmin:90:lonmax);
yticks(latmin:45:latmax);
end
cb = colorbar;
if handles.radiobutton3.Value == 1
   set(get(cb,'title'),'string','m^2/s^2');
elseif handles.radiobutton4.Value == 1
    set(get(cb,'title'),'string','mGal');
elseif handles.radiobutton5.Value == 1 
    set(get(cb,'title'),'string','E');
end
colormap(jet)

xlabel('Longitude(degree)');
ylabel('Latitude(degree)');

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,c]=uiputfile('*.txt','保存');
output=get(handles.uitable1,'data');
if c==1
    file=[pathname,filename];
    csvwrite(file,output);
    h=dialog('Name','Save data','Position',[200 200 200 70]);
    uicontrol('Style','text','Units','pixels','Position',[50 40 120 20],'FontSize',10,...
        'Parent',h,'String','save done');  
    uicontrol('Units','pixels','Position',[80 10 50 20],'FontSize',10,...
        'Parent',h,'String','OK','Callback','delete(gcf)');

end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable1,'data','');
set(handles.uitable2,'data','');

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_f_handle=figure('visible','off');
new_axes=copyobj(handles.axes1,new_f_handle); %picture是GUI界面绘图的坐标系句柄
set(new_axes,'units','default','Position','default');
 colorbar;
colormap(jet)
[filename,pathname]=uiputfile({'*.png'},'save picture as');
if ~filename
    return
else
    file=strcat(pathname,filename);
    print(new_f_handle,'-djpeg',file);
end
delete(new_f_handle);
h=dialog('Name','Save data','Position',[200 200 200 70]);
uicontrol('Style','text','Units','pixels','Position',[50 40 120 20],'FontSize',10,...
    'Parent',h,'String','save done'); 
uicontrol('Units','pixels','Position',[80 10 50 20],'FontSize',10,...
    'Parent',h,'String','OK','Callback','delete(gcf)');


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
cla reset
axis off


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('*.txt','open file');
if filename==0
    return
else
    file=[pathname,filename];
    data=csvread(file);
    datat=data(:,[1,2]);
set(handles.uitable7,'data',datat);

end 

% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
set(handles.edit1,'enable','on');
set(handles.edit2,'enable','on');
set(handles.edit3,'enable','on');
set(handles.edit4,'enable','on');
set(handles.edit5,'enable','on');
set(handles.edit6,'enable','on');
set(handles.pushbutton14,'enable','of');


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7
set(handles.edit1,'enable','of');
set(handles.edit2,'enable','of');
set(handles.edit3,'enable','of');
set(handles.edit4,'enable','of');
set(handles.edit5,'enable','of');
set(handles.edit6,'enable','of');
set(handles.pushbutton14,'enable','on');
