function [OTFParam,disp_map] = matrixstoring_2DG(xcenter,ycenter,xrange,yrange,imgfile,m)

% interchange X and Y to match with C++. X is row and Y is column
Xlb = round(xcenter)-xrange; % left boundary of x 
Xrb = round(xcenter)+xrange; % right boundary of x 
Ylob = round(ycenter)+yrange; % lower boundary of y 
Yub = round(ycenter)-yrange; % upper boundary of y 
[ysize,xsize]=size(imgfile);

% Check if the search range is within the range of the image
if Xrb>xsize
    disp('Center is near the bottom of the image!');
%     fprintf('The value of Xrb is %i at I_sumframe.cam%i(%i,%i,%i)\n',Xrb,m,i,j,k);
    Xrb=xsize;
    fprintf('The value of Xrb is set to %i\n',xsize)
end

if Xlb<1
    disp('Center is near the top of the image!');
%     fprintf('The value of Xlb is %i at I_sumframe.cam%i(%i,%i,%i)\n',Xlb,m,i,j,k);
    Xlb=1;
    fprintf('The value of Xlb is set to %i\n',1)
end
if Ylob>ysize
    disp('Center is near the right edge of the image!');
%     fprintf('The value of Ylob is %i at I_sumframe.cam%i(%i,%i,%i)\n',Ylob,m,i,j,k);
    Ylob=ysize;
    fprintf('The value of Ylob is set to %i\n',ysize)
end
if Yub<1
    disp('Center is near the left edge of the image!');
%     fprintf('The value of Yub is %i at I_sumframe.cam%i(%i,%i,%i)\n',Yub,m,i,j,k);
    Yub=1;
    fprintf('The value of Yub is set to %i\n',1)
end

x_mesh= Xlb:Xrb;
y_mesh = Yub:Ylob;
% [X,Y] = meshgrid(x_mesh,y_mesh);
zdata = double(imgfile(y_mesh,x_mesh)); % x and y are being interchanged to match with C++
OTFParam= fmgaussfit(x_mesh,y_mesh,zdata);

%% VSC disparities calculation 
% Equation: di=(dix,diy)=(xi',yi')-(xi,yi) ...
...or Mi'(X,Y,Z) = Mi(X,Y,Z)-di(X,Y,Z) where (xi',yi')is the center...
...of the particle from mapping functions and (xi,yi) is obtained from ...
...2D gaussian fitting %
dix=xcenter-OTFParam(5);
diy=ycenter-OTFParam(6);
disp_map=[dix,diy];



%% Storing matrix around the particles(not needed)
% nrow=Ylob-Yub+1;  % number of rows 
% ncol=Xrb-Xlb+1; % number of columns
% A=zeros(nrow,ncol);
% 
% for j = Yub:Ylob 
%     for i = Xlb:Xrb
%         A(i-Xlb+1,j-Yub+1)=tempfile(i,j);  
%     end  
% end

