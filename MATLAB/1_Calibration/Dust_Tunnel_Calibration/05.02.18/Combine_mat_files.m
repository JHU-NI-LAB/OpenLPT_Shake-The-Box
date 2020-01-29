load cam1.mat
camera_parameters(1)=camParaCalib;
clear camParaCalib
load cam2trans.mat
camera_parameters(2)=camParaCalib;
clear camParaCalib
load cam3trans.mat
camera_parameters(3)=camParaCalib;
clear camParaCalib
load cam4.mat
camera_parameters(4)=camParaCalib;
clear camParaCalib
load cam5.mat
camera_parameters(5)=camParaCalib;
clear camParaCalib
load cam6trans.mat
camera_parameters(6)=camParaCalib;
clear camParaCalib

camParaCalib=camera_parameters;
save('camParaCalib.mat','camParaCalib')