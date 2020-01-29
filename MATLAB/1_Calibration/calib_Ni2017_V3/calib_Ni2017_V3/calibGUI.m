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
    global img stop areanum centall lineall aa isx basepix baseidx gridspace

    [FileName,PathName] = uigetfile('*.*','Select the picture from paper');
    img=imread([PathName,FileName]);
    axes(handles.axes1);
    imshow(img);
    stop=0;
    areanum=0;
    centall=[];
    lineall=[];
    basepix=[];
    baseidx=[];
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
minline=str2num(mm{1})

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
    linecent  = flipud(linecent);
    linecent(:,3)=repmat(xpos,[size(linecent,1),1]);
    linecent(:,4)=repmat(0,[size(linecent,1),1]);
    
    totalnumlines=maxline-minline+1;

    
    
    if totalnumlines>size(linecent,1)
        incrematrix=sort(-diff(linecent(:,2)),'descend');
        meanincre=-mean(-incrematrix(totalnumlines-size(linecent,1)+1:end)); 
    
    end
    
    
    
    while totalnumlines>size(linecent,1)
        incre=-diff(linecent(:,2));
        fill=find(incre==max(incre));
        numfill=int8(max(incre)/meanincre)-1;
        numfill=double(numfill);
        insert=zeros(numfill+2,4);
        insert(:,3)=repmat(xpos,[size(insert,1),1]);
        insert(:,4)=repmat(0,[size(insert,1),1]);
        

        insert(:,1)=interp1([0,1+numfill],linecent(fill:fill+1,1),0:1+numfill);
        insert(:,2)=interp1([0,1+numfill],linecent(fill:fill+1,2),0:1+numfill);
        linecent = [linecent(1:fill-1,:);insert;linecent(fill+2:end,:)];
       
    end
    
    linecent(:,4)=minline:maxline;
    for ii=1:size(linecent,1)
         str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
         text(linecent(ii,1), linecent(ii,2), str, 'Color', 'g');hold on;
    end    
%     
%     
%     increment=sqrt(sum((linecent(end,:)-linecent(1,:)).^2))./(maxline-minline);
%     target=repmat(linecent(1,:),[size(linecent,1),1]);
%     dist=sqrt(sum((linecent(:,:)-target).^2,2));
%     index=maxline-round(dist./increment);
%     linecent(:,4)=index;
%     for ii=1:size(linecent,1)
%          str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
%          text(linecent(ii,1), linecent(ii,2), str, 'Color', 'k');hold on;
%     end
                   
else
%     linecent=sortrows(linecent,1);
%     linecent(:,4)=repmat(ypos,[size(linecent,1),1]);
%     increment=sqrt(sum((linecent(end,:)-linecent(1,:)).^2))./(maxline-minline);
%     target=repmat(linecent(1,:),[size(linecent,1),1]);
%     dist=sqrt(sum((linecent(:,:)-target).^2,2));
%     index=round(dist./increment)+minline;
%     linecent(:,3)=index;
%     for ii=1:size(linecent,1)
%          str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
%          text(linecent(ii,1), linecent(ii,2), str, 'Color', 'k');hold on;
%     end    
%     
    linecent=sortrows(linecent,1);
    linecent(:,4)=repmat(ypos,[size(linecent,1),1]);
    totalnumlines=maxline-minline+1;
    
    
    if totalnumlines>size(linecent,1)
        incrematrix=sort(diff(linecent(:,1)),'descend');
        meanincre=mean(incrematrix(totalnumlines-size(linecent,1)+1:end)); 
    end
    
    
    while totalnumlines>size(linecent,1)
        incre=diff(linecent(:,1));
        fill=find(incre==max(incre));
        numfill=int8(max(incre)/meanincre)-1;
        numfill=double(numfill);
        insert=zeros(numfill+2,4);
        insert(:,4)=repmat(ypos,[size(insert,1),1]);
        

        insert(:,1)=interp1([0,1+numfill],linecent(fill:fill+1,1),0:1+numfill);
        insert(:,2)=interp1([0,1+numfill],linecent(fill:fill+1,2),0:1+numfill);
        linecent = [linecent(1:fill-1,:);insert;linecent(fill+2:end,:)];
        
        
    end
    x_dir=input('What is the direction of +ve x? enter 1 for right or 2 for left')
    switch x_dir
        case 1
            linecent(:,3)=minline:maxline;
        case 2
            linecent(:,3)=maxline:-1:minline;
    end
            
%     if linecent(1,1)>linecent(size(linecent,1),1)
%     linecent(:,3)=minline:maxline;
%     else
%     linecent(:,3)=maxline:-1:minline;
%     end
   
    for ii=1:size(linecent,1)
         str = strcat('(',num2str(linecent(ii,3)),',',num2str(linecent(ii,4)),')');
         text(linecent(ii,1), linecent(ii,2), str, 'Color', 'g');hold on;
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
global finalnew
[fname,pth] = uiputfile('.txt'); % Type in name of file.
dlmwrite([pth,fname],finalnew,'precision','%.4f')


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global p2d3d holdon offseton gridspace noiterate

    [FileName,PathName] = uigetfile('*.*','4 cols or 5 cols');
    p2d3d=load([PathName,FileName]);
    if size(p2d3d,2)<5
    p2d3d=[p2d3d(:,1:2) zeros(size(p2d3d,1),1) p2d3d(:,3:4)];
    end
    holdon=0;
    offseton=0


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
global Npixw 
mm=get(hObject,'String');
Npixw=str2num(mm);


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
Npixh=str2num(mm);

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
hpix=str2num(mm);

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
wpix=str2num(mm);

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


%% --- calibration button
function pushbutton12_Callback(hObject, eventdata, handles)


% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Npixh Npixw hpix wpix p2d3d camParaCalib camParaCalib1 bb offseton gridspace noiterate

camParaknown.Npixh = Npixh;
camParaknown.Npixw = Npixw;
camParaknown.hpix = hpix;    % pixel size (mm)   [basler A504k = .012, mikrotron mc1631 = .014]
camParaknown.wpix = wpix;
set(handles.text10,'string','wait...');

%gridspace = handles.data;

p2d3d(:,3:5) = p2d3d(:,3:5).*gridspace;

    
    
camParaCalib1= db_calib_Tsai_Co_planar(p2d3d(:,1:2), p2d3d(:,3:5), camParaknown);
angle1 = rotm2eul(camParaCalib1.R)

if noiterate ==0
    if offseton==1
    [camParaCalib, angles] = iterate_calib(p2d3d(:,1:2), p2d3d(:,3:5), camParaCalib1);
    else 
    [camParaCalib, angles] = iterate_0offset(p2d3d(:,1:2), p2d3d(:,3:5), camParaCalib1);    
    end
else
    camParaCalib = camParaCalib1;
    angles = angle1;
end

set(handles.text10,'string',['least square fit result:' char(10) ...
    'f='   num2str(camParaCalib1.f_eff) ' ' ... 
    'offx=' num2str(camParaCalib1.Noffw) ' ' ... 
    'offy=' num2str(camParaCalib1.Noffh) ' ' ... 
    'errx=' num2str(camParaCalib1.err_x) ' '...
    'erry=' num2str(camParaCalib1.err_y) ...
    char(10) 'k1=' num2str(camParaCalib1.k1) ' ' ...
    'Euler Angles:' num2str(angle1(1)./pi.*180) ', ' num2str(angle1(2)./pi.*180) ', ' num2str(angle1(3)./pi.*180) ...
    char(10) char(10) ...
    'After nonlinear optimization:' char(10) ...
    'f='   num2str(camParaCalib.f_eff) ' ' ... 
    'offx=' num2str(camParaCalib.Noffw) ' ' ... 
    'offy=' num2str(camParaCalib.Noffh) ' ' ... 
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
global p2d3d camParaCalib1 camParaCalib holdon img
axes=handles.axes1;
if holdon==0
hold off;
end
imshow(img,[0 2000]); hold on;
% % plot camera view in 3D
% plot3(p2d3d(:,3).*20,p2d3d(:,4).*20,p2d3d(:,5).*20,'o');
% hold all;
% 
% pos1=[0,0,0];
% pos2=[-1,0,0];
% 
%     posnew1=pos1*camParaCalib.R'+camParaCalib.T';
%     posnew2=pos2*camParaCalib.R'+camParaCalib.T';
%     
%     quiver3(posnew1(1),posnew1(2),posnew1(3),posnew2(1)-posnew1(1),posnew2(2)-posnew1(2),posnew2(3)-posnew1(3),20);
%     
%     
%     
% pos1=[0,0,0];
% pos2=[-1,0,0];
% 
%     posnew1=pos1*camParaCalib1.R'+camParaCalib1.T';
%     posnew2=pos2*camParaCalib1.R'+camParaCalib1.T';
%     
%     plot3(posnew1(1),posnew1(2),posnew1(3),'r','MarkerFacecolor','r');
%     quiver3(posnew1(1),posnew1(2),posnew1(3),posnew2(1)-posnew1(1),posnew2(2)-posnew1(2),posnew2(3)-posnew1(3),20,'r');
%     
% set(handles.text10,'string',['arrows (blue arrow is the one after nonlinear iterations) indicate position and orientation of camera and circles indicate calibration target dots magnified by 20 times']);
%     


%% 

plot(p2d3d(:,1),p2d3d(:,2),'ro');hold on;

p2x=calibProj_Tsai(camParaCalib, p2d3d(:,3:5));
plot(p2x(:,1),p2x(:,2),'b+');

plot(p2d3d(:,1),p2d3d(:,2),'ro');
p2x=calibProj_Tsai(camParaCalib1, p2d3d(:,3:5));
plot(p2x(:,1),p2x(:,2),'gx');

set(handles.text10,'string',['red circle (image center), blue cross (after nonlinear optimization), green (least square fit)']);


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
global offseton
offseton=get(hObject,'Value');


% % --- Executes on button press in pushbutton16.
% function pushbutton16_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton16 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
%     global basepix baseidx bb centall xind yind
% 
%     clickxx=ginput(1);
%     if isempty(basepix)
%         basepix=clickxx;
%         baseidx=[xind,yind];
%     else
%         basepix=[basepix;clickxx];
%         tmp=[xind,yind];
%         baseidx=[baseidx;tmp];
%     end
%     
%     xc= centall(:,1);
%     yc= centall(:,2);
%     
%     dist = (xc-clickxx(1)).^2+(yc-clickxx(2)).^2;
%     [mindist i0] = min(dist);
% 
%     plot(xc(i0), yc(i0), 'rx','markers',12);
%     hold on;
% 
% 
% function edit11_Callback(hObject, eventdata, handles)
% % hObject    handle to edit11 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit11 as text
% %        str2double(get(hObject,'String')) returns contents of edit11 as a double
% global xind 
% mm=get(hObject,'String');
% xind=str2num(mm{1});
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit11_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit11 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% 
% function edit12_Callback(hObject, eventdata, handles)
% % hObject    handle to edit12 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit12 as text
% %        str2double(get(hObject,'String')) returns contents of edit12 as a double
% global yind 
% mm=get(hObject,'String');
% yind=str2num(mm{1});
% 
% % --- Executes during object creation, after setting all properties.
% function edit12_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit12 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% % --- Executes on button press in pushbutton17.
% function pushbutton17_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton17 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
%     global basepix baseidx bb centall xind yind centall
% 
%     x0 = basepix(1,1);
%     y0 = basepix(1,2);
%     i0xind = baseidx(1,1);
%     i0yind = baseidx(1,2);
%     
% 
%     x1 = basepix(2,1);
%     y1 = basepix(2,2);
%     i1xind = baseidx(2,1);
%     i1yind = baseidx(2,2);
%         
%     x2 = basepix(3,1);
%     y2 = basepix(3,2);
%     i2xind = baseidx(3,1);
%     i2yind = baseidx(3,2);
%     
%     xc = centall(:,1);
%     yc = centall(:,2);
% 
%     Np = size(xc,1);
%     
%     dist = (xc-x0).^2+(yc-y0).^2;
%     [mindist i0] = min(dist);
% 
%     dist = (xc-x1).^2+(yc-y1).^2;
%     [mindist i1] = min(dist);
% 
%     dist = (xc-x2).^2+(yc-y2).^2;
%     [mindist i2] = min(dist);
%         
% 
%     % Need to know a point close to the camera in order to
%     % compensate for the perspective projection distortion when find particle coordinates
%     set(handles.text10,'string',['click a point that is the closest to the camera']);
%     but = 0;
%     while but ~= 1
%         [xcam0 ycam0 but] = ginput(1);
%     end
% 
%         
%         % Now determine the point coordinates
%         % first, form two base vectors on the mask
%         e1 = [i1xind-i0xind, i1yind-i0yind];
%         e2 = [i2xind-i0xind, i2yind-i0yind];
%         % The prjection of these two vectors on image plane
%         e1p = [xc(i1)-xc(i0), yc(i1)-yc(i0)];
%         e2p = [xc(i2)-xc(i0), yc(i2)-yc(i0)];
%         e1pnorm = sum(e1p.^2);
%         e2pnorm = sum(e2p.^2);
%         e1pe2p = sum(e1p.*e2p);
%         d = (e1pnorm*e2pnorm - e1pe2p*e1pe2p);		% The denominator
%         % calaulte the coords of all points using the two base vectors
%         pind = zeros(Np, 2);
%         lenref = (xcam0-xc(i0))^2+(ycam0-yc(i0))^2;
%         for i=1:Np
%             c = [xc(i)-xc(i0), yc(i)-yc(i0)];
%             % The coefficient is used to compensate the effect of perspective projection
%             coef = 1 + 0.0*(((xc(i)-xcam0)^2+(yc(i)-ycam0)^2) - sum(c.^2))/lenref;
%             A = coef*(sum(c.*e1p)*e2pnorm - sum(c.*e2p)*e1pe2p)/d;
%             B = coef*(sum(c.*e2p)*e1pnorm - sum(c.*e1p)*e1pe2p)/d;
%             pind(i,:)=[A B];
%         end
%         % Now calculate the two components of dots' 3D coordinates on the mask plane
%         pmask = zeros(Np,2);
%         for i=1:Np
%             pmask(i,:) = floor(e1*pind(i,1) +0.5) + floor(e2*pind(i,2)+0.5) + [i0xind i0yind];
%         end
%         for i=1:Np
%             str = strcat('(',num2str(pmask(i,1)),',',num2str(pmask(i,2)),')');
%             text(xc(i), yc(i), str, 'Color', 'r');
%         end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global isx xpos ypos xx lineall maxline minline centall finalnew isextra noiterate
% para matrix col 1 slop col 2 interception col3 indx col4 indy
% fit a linear curve over x and y directions

ll=sortrows(lineall,4);
[c,ia,ic]=unique(ll(:,4),'legacy');
iidx=diff(ia);
twolines=find(iidx==2);
twolines=ia(twolines(1:end))+1;
k=1;


for lii=1:length(twolines)

    idxstep=sign(ll(twolines(lii)+1,3)-ll(twolines(lii),3));

    parax(k,4)=ll(twolines(lii),4);
    parax(k,3)=0;
    
    parax(k,1:2)=polyfit([ll(twolines(lii),1),ll(twolines(lii)+1,1)],[ll(twolines(lii),2),ll(twolines(lii)+1,2)],1);
%     
%     tmpxx=zeros(length((ll(twolines(lii),3):idxstep:ll(twolines(lii)+1,3))),4);
%     tmpxx(:,3)=(ll(twolines(lii),3):idxstep:ll(twolines(lii)+1,3))';
%     step=idxstep.*(ll(twolines(lii)+1,1)-ll(twolines(lii),1))./(ll(twolines(lii)+1,3)-ll(twolines(lii),3));
%     tmpxx(:,1)=(ll(twolines(lii),1):step:ll(twolines(lii)+1,1))';
% 
%     step=idxstep.*(ll(twolines(lii)+1,2)-ll(twolines(lii),2))./(ll(twolines(lii)+1,3)-ll(twolines(lii),3));
%     tmpxx(:,2)=(ll(twolines(lii),2):step:ll(twolines(lii)+1,2))';
%     tmpxx(:,4)=ones(size(tmpxx,1),1).*ll(twolines(lii),4);
%     if lii==1
%         allnew=tmpxx;
%     else 
%         allnew=[allnew;tmpxx];
%     end
    k=k+1;
end




ll=sortrows(lineall,3);
[c,ia,ic]=unique(ll(:,3),'legacy');
iidx=diff(ia);
twolines=find(iidx==2);
twolines=ia(twolines(1:end))+1;

k=1;
for lii=1:length(twolines)

    idxstep=sign(ll(twolines(lii)+1,4)-ll(twolines(lii),4));

    paray(k,3)=ll(twolines(lii),3);
    paray(k,4)=0;
   
    paray(k,1:2)=polyfit([ll(twolines(lii),1),ll(twolines(lii)+1,1)],[ll(twolines(lii),2),ll(twolines(lii)+1,2)],1);
    k=k+1;

end

% 
% 
% parax_new(:,4)=min(parax(:,4)).*2:max(parax(:,4)).*2;
% tmp1=polyfit(parax(:,4),parax(:,1),1);
% parax_new(:,1)=tmp1(1).*parax_new(:,4)+tmp1(2);
% 
% tmp1=polyfit(parax(:,4),parax(:,2),1);
% parax_new(:,2)=tmp1(1).*parax_new(:,4)+tmp1(2);
% 
% 
% paray_new(:,3)=min(paray(:,3)).*2:max(paray(:,3)).*2;
% tmp1=polyfit(paray(:,3),paray(:,1),1);
% paray_new(:,1)=tmp1(1).*paray_new(:,3)+tmp1(2);
% 
% tmp1=polyfit(paray(:,3),paray(:,2),1);
% paray_new(:,2)=tmp1(1).*paray_new(:,3)+tmp1(2);

parax_new(:,4)=min(parax(:,4)).*2:max(parax(:,4)).*2;
paray_new(:,3)=min(paray(:,3)).*2:max(paray(:,3)).*2;


if isextra==1
    parax_new(:,1)=interp1(parax(:,4),parax(:,1),min(parax(:,4)).*2:max(parax(:,4)).*2,'linear','extrap');
    parax_new(:,2)=interp1(parax(:,4),parax(:,2),min(parax(:,4)).*2:max(parax(:,4)).*2,'linear','extrap');

    paray_new(:,1)=interp1(paray(:,3),paray(:,1),min(paray(:,3)).*2:max(paray(:,3)).*2,'linear','extrap');
    paray_new(:,2)=interp1(paray(:,3),paray(:,2),min(paray(:,3)).*2:max(paray(:,3)).*2,'linear','extrap');
else
    parax_new=parax;
    paray_new=paray;
end


k=1
for i=1:size(parax_new,1)
    for j=1:size(paray_new,1)
        allnew(k,3)=paray_new(j,3);
        allnew(k,4)=parax_new(i,4);
        pixx=(paray_new(j,2)-parax_new(i,2))./(parax_new(i,1)-paray_new(j,1));
        pixy=paray_new(j,1).*pixx+paray_new(j,2);
        allnew(k,1)=pixx;
        allnew(k,2)=pixy;
        k=k+1;
    end
end





searchrad=10;
k=1;
finalnew=[];
for inew=1:size(allnew,1)
    xc = allnew(inew,1);
    yc = allnew(inew,2);
    dist = sqrt((xc-centall(:,1)).^2+(yc-centall(:,2)).^2);
    [mindist i0] = min(dist);
    if mindist<searchrad
        finalnew(k,:)=[centall(i0,:),allnew(inew,3:4)];
        k=k+1;
        str = strcat('(',num2str(allnew(inew,3)),',',num2str(allnew(inew,4)),')');
        text(centall(i0,1), centall(i0,2), str, 'Color', 'r');hold on;
    
    end
end

plot(allnew(:,1),allnew(:,2),'rx');
plot(finalnew(:,1),finalnew(:,2),'m*');


    



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

global gridspace
mm=get(hObject,'String');
gridspace=str2num(mm);
% handles.data=gridspace;
% guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
global isextra
isextra=get(hObject,'Value');


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
global noiterate
noiterate=get(hObject,'Value');