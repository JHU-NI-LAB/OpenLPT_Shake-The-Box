dir='./05.02.18/';
cam=6;
centdistx=43;
centdisty=10;
gridspace_this=6;
gridspace_that=3;
%%
a=load([dir 'cam' num2str(cam,'%01i') '.txt']);
a(:,3)=(centdistx+(a(:,3).*gridspace_this))./gridspace_that;
a(:,4)=(centdisty+(a(:,4).*gridspace_this))./gridspace_that;
dlmwrite([dir 'cam' num2str(cam,'%01i') 'trans.txt'],a,'precision','%.4f');