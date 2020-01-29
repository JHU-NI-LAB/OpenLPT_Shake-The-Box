function dispar_brightest = VSCpeakfinder(exsize,eysize,stepsize,VSC_xrange,VSC_yrange,Isub_sumframe,m,i,j,k)
% dx=-epsln_x:stepsize:epsln_x;
% dy=-epsln_y:stepsize:epsln_y;
zData = double(cell2mat(Isub_sumframe));
% if there is no particle in that subvolume, disparity is not available (NA) 
if max(zData(:))== 0
    dispar_brightest = ['NA,', 'NA'];
    fprintf('No particle at I_sumframe.cam%i(%i,%i,%i)\n',m,i,j,k)
    return
end
[~, ind] = max(zData(:)); % amp is the amplitude.
[xcenter,ycenter] = ind2sub(size(zData),ind);

Xlb = round(xcenter)-VSC_xrange;  % left boundary of x
Xrb = round(xcenter)+ VSC_xrange; % right boundary of x 
Ylob = round(ycenter)+ VSC_yrange;     % lower boundary of y 
Yub = round(ycenter)-VSC_yrange;      % upper boundary of y 
[xsize,ysize]=size(zData);

% Check if the search range is within the range of the image
if Xrb>xsize
    disp('Center is near the bottom of the image!');
    fprintf('The value of Xrb is %i at I_sumframe.cam%i(%i,%i,%i)\n',Xrb,m,i,j,k);
    Xrb=xsize;
    fprintf('The value of Xrb is set to %i\n',xsize)
end
if Xlb<1
    disp('Center is near the top of the image!');
    fprintf('The value of Xlb is %i at I_sumframe.cam%i(%i,%i,%i)\n',Xlb,m,i,j,k);
    Xlb=1;
    fprintf('The value of Xlb is set to %i\n',1)
end
if Ylob>ysize
    disp('Center is near the right edge of the image!');
    fprintf('The value of Ylob is %i at I_sumframe.cam%i(%i,%i,%i)\n',Ylob,m,i,j,k);
    Ylob=ysize;
    fprintf('The value of Ylob is set to %i\n',ysize)
end
if Yub<1
    disp('Center is near the left edge of the image!');
    fprintf('The value of Yub is %i at I_sumframe.cam%i(%i,%i,%i)\n',Yub,m,i,j,k);
    Yub=1;
    fprintf('The value of Yub is set to %i\n',1)
end
x_mesh= Xlb:Xrb;
y_mesh = Yub:Ylob;
z_gaussian = zData(x_mesh,y_mesh); % Does not need to interchange x and y because the index is from matlab 
VSCresults = fmgaussfit(x_mesh,y_mesh,z_gaussian);
xc=(VSCresults(5)-(exsize/2+0.5))*stepsize; % Convert matlab disparity position to actual disparty value
yc=(VSCresults(6)-(eysize/2+0.5))*stepsize; % Convert matlab disparity position to actual disparty value

dispar_brightest = [xc, yc];