% load camParaCalib4
function CameraDistribution(camParaCalib)
camera_parameters=camParaCalib(1:4);
wo=[0 0 0]';
wx=[100 0 0]';
wy=[0 100 0]';
wz=[0 0 100]';
co=[0 0 0]'; %camera origin in camera coordinate
cdx=[100 0 0]';
cdy=[0 100 0]';
cdz=[0 0 100]';
% figure;
for cam=1:4

plot3([wo(1);wx(1)],[wo(2);wx(2)],[wo(3);wx(3)],'b','LineWidth',4);hold on
plot3([wo(1);wy(1)],[wo(2);wy(2)],[wo(3);wy(3)],'k','LineWidth',4);hold on
plot3([wo(1);wz(1)],[wo(2);wz(2)],[wo(3);wz(3)],'r','LineWidth',4);hold on

cow=camera_parameters(cam).R\(co-camera_parameters(cam).T);
plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'k');hold on

cdxw=camera_parameters(cam).R\(cdx-camera_parameters(cam).T);
plot3([cow(1);cdxw(1)],[cow(2);cdxw(2)],[cow(3);cdxw(3)],'b','LineWidth',4);hold on

cdyw=camera_parameters(cam).R\(cdy-camera_parameters(cam).T);
plot3([cow(1);cdyw(1)],[cow(2);cdyw(2)],[cow(3);cdyw(3)],'k','LineWidth',4);hold on

cdzw=camera_parameters(cam).R\(cdz-camera_parameters(cam).T);
plot3([cow(1);cdzw(1)],[cow(2);cdzw(2)],[cow(3);cdzw(3)],'r','LineWidth',4);hold on

scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
scatter3([0],[0],[0],'r*')
axis equal
end
end

% cow=camParaCalib(2).R\(co-camParaCalib(2).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camParaCalib(2).R\(cdz-camParaCalib(2).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
% 
% cow=camParaCalib(3).R\(co-camParaCalib(3).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camParaCalib(3).R\(cdz-camParaCalib(3).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
% 
% cow=camParaCalib(4).R\(co-camParaCalib(4).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camParaCalib(4).R\(cdz-camParaCalib(4).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
% 
% cow=camera_parameters(5).R\(co-camera_parameters(5).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camera_parameters(5).R\(cdz-camera_parameters(5).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
% 
% cow=camera_parameters(6).R\(co-camera_parameters(6).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camera_parameters(6).R\(cdz-camera_parameters(6).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')
% 
% cow=camera_parameters(7).R\(co-camera_parameters(7).T);
% plot3([cow(1);wo(1)],[cow(2);wo(2)],[cow(3);wo(3)],'b');hold on
% cdw=camera_parameters(7).R\(cdz-camera_parameters(7).T);
% plot3([cow(1);cdw(1)],[cow(2);cdw(2)],[cow(3);cdw(3)],'r','LineWidth',4);hold on
% scatter3([cow(1)],[cow(2)],[cow(3)],'g*')


% x=[0 0 0]'
% X=[camParaCalib(1).R*x+camParaCalib(1).T]
% plot3([x(1);X(1)],[x(2);X(2)],[x(3);X(3)],'k');hold on
% % x=camParaCalib(1).R\(X-camParaCalib(1).T)
% Z=[0 0 1]'
% z=camParaCalib(1).R\(Z-camParaCalib(1).T)
% plot3([X(1);-Z(1)],[X(2);-Z(2)],[X(3);-Z(3)],'b');hold on
% 
% X=-[camParaCalib(2).R*x+camParaCalib(2).T];
% plot3([x(1);X(1)],[x(2);X(2)],[x(3);X(3)],'k');hold on
% X=-[camParaCalib(3).R*x+camParaCalib(3).T];
% plot3([x(1);X(1)],[x(2);X(2)],[x(3);X(3)],'k');hold on
% X=-[camParaCalib(4).R*x+camParaCalib(4).T];
% plot3([x(1);X(1)],[x(2);X(2)],[x(3);X(3)],'k');hold on

% Xtest3D=[0 0 0];
% Xc = Xtest3D * (camParaCalib(1).R)';
% Xc(:,1) = Xc(:,1) + camParaCalib(1).T(1);
% Xc(:,2) = Xc(:,2) + camParaCalib(1).T(2);
% Xc(:,3) = Xc(:,3) + camParaCalib(1).T(3);
% plot3([x(1);X(1)],[x(2);X(2)],[x(3);X(3)],'k');hold on