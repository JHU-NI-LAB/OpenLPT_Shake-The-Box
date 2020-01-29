function [camParaCalib_new angles]=iterate_calib(p2d,p3d,camParaCalib)
%%%Calibration

camParaCalib_new=camParaCalib;
camParaCalib_new.kstar=1;

%% fine tune %%
[angles rotout] = rni_rotmat2angles(camParaCalib.R)

%params = [angles camParaCalib.T' camParaCalib.f_eff camParaCalib.k1 camParaCalib.Noffw camParaCalib.Noffh];
params = [angles camParaCalib.T' camParaCalib.f_eff camParaCalib.k1];
params_const = [camParaCalib.wpix, camParaCalib.hpix, camParaCalib.Npixw, camParaCalib.Npixh];
    
x=params;
options = optimset('MaxFunEvals',1e7,'MaxIter',1e7);
[x,fval,exitflag,output] = fminsearch( @ftsai_min_0offset, params, options, params_const, p3d, p2d)

for i=1:5
[x,fval,exitflag,output] = fminsearch( @ftsai_min_0offset, x, options, params_const, p3d, p2d)
end
%[x,fval,exitflag,output] = fminsearch( @ftsai_min, x, options, params_const, pall(:,1:3), pall(:,4:5))


%[x,fval,exitflag,output] = patternsearch( @(x)ftsai_min(x, params_const, pall(:,1:3), pall(:,4:5)), x, [])

%%
%[x,fval,exitflag,output] = ga( @(x)ftsai_min(x, params_const, pall(:,1:3), pall(:,4:5)), length(params), [],[],[],[],params-0.1.*abs(params),params+0.1.*abs(params))

params=x;

    theta=params(1);
    phi=params(2);
    mu=params(3);
    angles=params(1:3);
    
    rotztheta=[cos(theta), -sin(theta), 0; +sin(theta), cos(theta), 0; 0 0 1];
    rotzmu=[cos(mu), -sin(mu), 0; +sin(mu), cos(mu), 0; 0 0 1];
    rotyphi=[cos(phi), 0, sin(phi); 0 1 0; -sin(phi), 0, cos(phi)];
    reflect = [1 0 0;0 1 0;0 0 -1];

    camParaCalib_new.R = reflect*rotztheta*rotyphi*rotzmu;

    camParaCalib_new.T = params(4:6)';
    
    
    
    T = [[camParaCalib_new.R params(4:6)']; [0 0 0 1]];
    Tinv = T^(-1);
    camParaCalib.Rinv = Tinv(1:3, 1:3);
    camParaCalib.Tinv = Tinv(1:3, 4);



    camParaCalib_new.f_eff = params(7);
    camParaCalib_new.k1 = params(8);
%     camParaCalib_new.Noffw = params(9);
%     camParaCalib_new.Noffh = params(10);
    camParaCalib_new.Noffw = 0;
    camParaCalib_new.Noffh = 0;
    
    
    camParaCalib_new.Rinv = Tinv(1:3, 1:3);
    camParaCalib_new.Tinv = Tinv(1:3, 4);
    

    
        Xproj = calibProj_Tsai(camParaCalib_new, p3d);
        err_img = Xproj - p2d;
        camParaCalib_new.err_x = sqrt(mean((err_img(:,1)).^2));
        camParaCalib_new.err_y = sqrt(mean((err_img(:,2)).^2));
        camParaCalib_new.err_t = sqrt((camParaCalib_new.err_x)^2 + (camParaCalib_new.err_y)^2);  
end