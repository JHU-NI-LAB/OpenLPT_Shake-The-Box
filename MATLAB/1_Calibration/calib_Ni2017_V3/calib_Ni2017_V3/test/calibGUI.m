function varargout = calibGUI(varargin)
% CALIBGUI MATLAB code for calibGUI.fig
%      CALIBGUI, by itself, creates a new CALIBGUI or raises the existing
%      singleton*.
%
%      H = CALIBGUI returns the handle to a new CALIBGUI or the handle to
%      the existing singleton*.
%
%      CALIBGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBGUI.M with the given input arguments.
%
%      CALIBGUI('Property','Value',...) creates a new CALIBGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibGUI

% Last Modified by Ni v1 14-Feb-2017 17:27:47
% The first version:
%
% 

% Begin initialization code - DO NOT EDIT

global img part_radius stop


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @calibGUI_OutputFcn, ...
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

% --- Executes just before calibGUI is made visible.
function calibGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibGUI (see VARARGIN)

% Choose default command line output for calibGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calibGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in Open button to open image file.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global img stop areanum centall lineall aa isx

    [FileName,PathName] = uigetfile('*.*','Select the picture from paper');
    img=imread([PathName,FileName]);
    axes(handles.axes1);
    imshow(img);
    stop=0;
    areanum=0;
    centall=[];
    lineall=[];
    aa=[];
    isx=1;
    set(handles.text10,'string','Fill in the dot rad (How many pixels of the radius for the smallest dot)');
    


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global img part_radius bb cent centall

cut=get(hObject,'Value');
bb=im2bw(img,cut./255);
axes(handles.axes1);
imshow(bb);

b1=imcomplement(bb);
s = regionprops(b1,'centroid','Area');
centroids = cat(1, s.Centroid);
areas = cat(1,s.Area);

hold on
cent=centroids(find(areas>part_radius));
plot(centroids(:,1),centroids(:,2), 'r+');

if ~isempty(centall)
    plot(centall(:,1),centall(:,2), 'bo');
end
set(handles.text10,'string','Now you can remove points by click remove points and use cursor to draw a rectangle to remove alll dots inside. You can repeat this for many times. When you are done, you can click Add area. If you have multiple areas, always remove all points for the previous areas shown in blue circles.');



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bb cent part_radius stop centall


axes(handles.axes1);

rect = getrect;
hold off;

imshow(bb);
minl=max(1,rect(2));
maxl=min(rect(2)+rect(4),size(bb,1));

minr=max(1,rect(1));
maxr=min(rect(1)+rect(3),size(bb,2));

bb(minl:maxl,minr:maxr)=0;

b1=imcomplement(bb);
s = regionprops(b1,'centroid','Area');
centroids = cat(1, s.Centroid);
areas = cat(1,s.Area);

hold on
cent=centroids(find(areas>part_radius),:);
plot(cent(:,1),cent(:,2), 'r+');

if ~isempty(centall)
    plot(centall(:,1),centall(:,2), 'bo');
end
set(handles.text10,'string','Remove wrong dataspoints by draw a rectangle using your cursor.');


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global part_radius;
num=get(hObject,'String');
part_radius=str2num(num);

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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global centall cent
if ~isempty(centall)
centall=[centall;cent];
else
    centall=cent;
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global centall bb
rect = getrect;

indx=find(centall(:,1)<rect(1)+rect(3)&centall(:,1)>rect(1)&centall(:,2)<rect(2)+rect(4)&centall(:,2)>rect(2));
centall(indx,:)=[];

imshow(bb)
plot(centall(:,1),centall(:,2), 'bo');


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
        global linenum xpos ypos isx aa
         contents = cellstr(get(hObject,'String'));
         aa=contents{get(hObject,'Value')};

         switch aa;
             
         case 'x' % User selects Peaks.
            xpos=linenum;
            isx=1;
            ypos=[];
            
         case 'y' % User selects Membrane.
            ypos=linenum;
            isx=0;
            xpos=[];
       
         end
         
         


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global linenum aa xpos isx ypos
mm=get(hObject,'String');
linenum=str2num(mm{1});

if ~isempty(aa)
         switch aa;
             
         case 'x' % User selects Peaks.
            xpos=linenum;
            isx=1;
            ypos=[];
            
         case 'y' % User selects Membrane.
            ypos=linenum;
            isx=0;
            xpos=[];
       
         end
else 
            xpos=linenum;
            isx=1;
            ypos=[];
end
         
         
         

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
global maxline 
mm=get(hObject,'String');
maxline=str2num(mm{1});

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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global minline 
mm=get(hObject,'String');
minline=str2num(mm{1});

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


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global linecent isx xpos ypos xx lineall maxline minline centall bb
dist=get(hObject,'Value');
hold off;
imshow(bb);hold on;
plot(centall(:,1),centall(:,2),'bo');

plot(xx(:,1),xx(:,2),'r');
distance = point_to_line(centall,xx(1,:),xx(2,:));
plot(centall(distance<dist,1),centall(distance<dist,2),'ro','MarkerFacecolor','r')

linecent=centall(find(distance<dist),:);

    
    

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xx bb centall

xx=ginput(2);

hold off;

imshow(bb); hold on; 
plot(centall(:,1),centall(:,2),'bo');
plot(xx(:,1),xx(:,2),'r');


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global linecent isx xpos ypos xx lineall maxline minline centall

if isx==1
    linecent=sortrows(linecent,2);
    linecent(:,3)=repmat(xpos,[size(linecent,1),1]);
    increment=sqrt(sum((linecent(end,:)-linecent(1,:)).^2))./(maxline-minline);
    target=repmat(linecent(1,:),[size(linecent,1),1]);
    dist=sqrt(sum((linecent(:,:)-target).^2,2));
    index=maxline-round(dist./increment);
    linecent(:,4)=index;
    for ii=1:size(linecent,1)
         str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
         text(linecent(ii,1), linecent(ii,2), str, 'Color', 'k');hold on;
    end
                   
else
    linecent=sortrows(linecent,1);
    linecent(:,4)=repmat(xpos,[size(linecent,1),1]);
    increment=sqrt(sum((linecent(end,:)-linecent(1,:)).^2))./(maxline-minline);
    target=repmat(linecent(1,:),[size(linecent,1),1]);
    dist=sqrt(sum((linecent(:,:)-target).^2,2));
    index=round(dist./increment)+minline;
    linecent(:,3)=index;
    for ii=1:size(linecent,1)
         str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
         text(linecent(ii,1), linecent(ii,2), str, 'Color', 'k');hold on;
    end    
    
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lineall linecent
if isempty(lineall)
    lineall=linecent;
else
    lineall=[lineall;linecent];
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lineall
[fname,pth] = uiputfile('.txt'); % Type in name of file.
dlmwrite([pth,fname],lineall,'precision','%.4f')


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global p2d3d holdon

    [FileName,PathName] = uigetfile('*.*','4 cols or 5 cols');
    p2d3d=load([PathName,FileName]);
    if size(p2d3d,2)<5
    p2d3d=[p2d3d(:,1:2) zeros(size(p2d3d,1),1) p2d3d(:,3:4)];
    end
    holdon=0;
    



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
global Npixw 
mm=get(hObject,'String');
Npixw=str2num(mm{1});


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global Npixh 
mm=get(hObject,'String');
Npixh=str2num(mm{1});

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
global hpix
mm=get(hObject,'String');
hpix=str2num(mm{1});

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
global wpix 
mm=get(hObject,'String');
wpix=str2num(mm{1});

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Npixh Npixw hpix wpix p2d3d camParaCalib camParaCalib1 bb
camParaknown.Npixh = Npixh;
camParaknown.Npixw = Npixw;
camParaknown.hpix = hpix;    % pixel size (mm)   [basler A504k = .012, mikrotron mc1631 = .014]
camParaknown.wpix = wpix;
set(handles.text10,'string','wait...');
camParaCalib1= db_calib_Tsai_Co_planar(p2d3d(:,1:2), p2d3d(:,3:5), camParaknown);
angle1 = rni_rotmat2angles(camParaCalib1.R);

[camParaCalib, angles] = iterate_calib(p2d3d(:,1:2), p2d3d(:,3:5), camParaCalib1);
set(handles.text10,'string',['least square fit result:' char(10) ...
    'f='   num2str(camParaCalib1.f_eff) char(10) ... 
    'errx=' num2str(camParaCalib1.err_x) ' '...
    'erry=' num2str(camParaCalib1.err_y) ...
    char(10) 'k1=' num2str(camParaCalib1.k1) ' ' ...
    'Euler Angles:' num2str(angle1(1)./pi.*180) ', ' num2str(angle1(2)./pi.*180) ', ' num2str(angle1(3)./pi.*180) ...
    char(10) char(10) ...
    'After nonlinear optimization:' char(10) ...
    'f='   num2str(camParaCalib.f_eff) char(10) ... 
    'errx=' num2str(camParaCalib.err_x) ' '...
    'erry=' num2str(camParaCalib.err_y) ' ' ...
    char(10) 'k1=' num2str(camParaCalib.k1) ' ' ...  
    'Euler Angles:' num2str(angles(1)./pi.*180) ', ' num2str(angles(2)./pi.*180) ', ' num2str(angles(3)./pi.*180)]);




% 
% 
% axes(handles.axes1);
% 
% hold off;
% imshow(bb);hold on;
% plot(p2d3d(:,1),p2d3d(:,2), 'bo');
% p2d = calibProj_Tsai(camParaCalib, p2d3d(:,3:5));
% plot(p2d(:,1),p2d(:,2), 'r+');
% 

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global camParaCalib
[fname,pth] = uiputfile('.mat'); % Type in name of file.
save([pth,fname],'camParaCalib');


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p2d3d camParaCalib1 camParaCalib holdon
axes=handles.axes1;
if holdon==0
hold off;
end

plot3(p2d3d(:,3).*20,p2d3d(:,4).*20,p2d3d(:,5).*20,'o');
hold all;

pos1=[0,0,0];
pos2=[-1,0,0];

    posnew1=pos1*camParaCalib.R'+camParaCalib.T';
    posnew2=pos2*camParaCalib.R'+camParaCalib.T';
    
    quiver3(posnew1(1),posnew1(2),posnew1(3),posnew2(1)-posnew1(1),posnew2(2)-posnew1(2),posnew2(3)-posnew1(3),20);
    
    
    
pos1=[0,0,0];
pos2=[-1,0,0];

    posnew1=pos1*camParaCalib1.R'+camParaCalib1.T';
    posnew2=pos2*camParaCalib1.R'+camParaCalib1.T';
    
    plot3(posnew1(1),posnew1(2),posnew1(3),'r','MarkerFacecolor','r');
    quiver3(posnew1(1),posnew1(2),posnew1(3),posnew2(1)-posnew1(1),posnew2(2)-posnew1(2),posnew2(3)-posnew1(3),20,'r');
    
set(handles.text10,'string',['arrows (blue arrow is the one after nonlinear iterations) indicate position and orientation of camera and circles indicate calibration target dots magnified by 20 times']);
    
    %posnew=posnew+camParaCalib.T;
    %posnew=posnew'
%     
% %     allcol=[posnew(:,1)+center(1),posnew(:,2)+center(2),posnew(:,3)+center(3)];
% %     all1=sortrows(allcol,[1,2]);
%     xx1=reshape(posnew(:,1)+camParaCalib.T(1),[2,(n+1)]);
%     yy1=reshape(posnew(:,2)+camParaCalib.T(2),[2,(n+1)]);
%     zz1=reshape(posnew(:,3)+camParaCalib.T(3),[2,(n+1)]);
% 
% surf(xx1,yy1,zz1.*20);
