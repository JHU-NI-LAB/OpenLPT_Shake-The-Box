function correct_calib(fileload, filesave, camParaCalib, gridspace, platethick)


%% load file with target dot 2D and 3D positions
cam=load (fileload);

cam3d = [zeros(size(cam,1),1),cam(:,3:4)].*gridspace;




%% switch center

cam2d=cam(:,1:2);
cam2d(:,1)=(cam2d(:,1)-camParaCalib.Npixw/2-camParaCalib.Noffw)*camParaCalib.wpix;
cam2d(:,2)=(-cam2d(:,2)+camParaCalib.Npixh/2-camParaCalib.Noffh)*camParaCalib.hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative


%% organize calibration file
calarray=mat2calarray(camParaCalib);
R(:,:) = eul2rotm(calarray(1:3));
Rinv(:,:) = inv(R(:,:));
for ii=1:size(cam2d,1)
    u(ii,:) = img2unitv(calarray(:), Rinv(:,:), cam2d(ii,:));
end

%% get direction

theta1=acos(u(:,1));
theta1(theta1>pi/2)=pi-theta1(theta1>pi/2);
theta2=asin(sin(theta1)*1.33./1.49)

u2 = cos(theta2);
vert = cross(u,repmat([1,0,0],[size(u,1) 1]));
vert = vert./repmat(sqrt(sum(vert.^2,2)),[1,3]);

cam3d_new=cam3d;

for ii=1:size(u,1)
    b = 1-u2(ii,1).^2;
    c = vert(ii,2);
    d = vert(ii,3);
    
    y1=sqrt(b/((d/c)^2+1));
    x1 = -y1*(d/c);
%     syms x1 y1
%     eqns = [x1^2 + y1^2 == b, x1*c + y1*d == 0];
%     S = solve(eqns, [x1 y1]);
    u2(ii,2:3)=[x1,y1];
    disp(ii) = platethick.*tan(theta1(ii))-platethick.*tan(theta2(ii));
    u2d=u(ii,2:3);
    cam3d_new(ii,2) = cam3d(ii,2)+u2d(1)./norm(u2d).*disp(ii);
    cam3d_new(ii,3) = cam3d(ii,3)+u2d(2)./norm(u2d).*disp(ii);
    
end

cam=load (fileload);
cam(:,3:4)=cam3d_new(:,2:3)./gridspace;
dlmwrite(filesave,cam,'precision','%.4f')
end