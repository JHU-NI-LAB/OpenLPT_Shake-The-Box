function [camParaCalib_new angles]=iterate_calib(p2d,p3d,camParaCalib)
    %%%Calibration

    camParaCalib_new=camParaCalib;
    camParaCalib_new.k1star=0;

    %% fine tune %%
    eul = rotm2eul(camParaCalib.R);
    angles = eul';
    
    params = [angles camParaCalib.T' camParaCalib.f_eff camParaCalib.k1 camParaCalib.Noffw camParaCalib.Noffh];
    params_const = [camParaCalib.wpix, camParaCalib.hpix, camParaCalib.Npixw, camParaCalib.Npixh];

    x=params;
    options = optimset('MaxFunEvals',1e7,'MaxIter',1e7);
    [x,fval,exitflag,output] = fminsearch( @ftsai_min, params, options, params_const, p3d, p2d)

    for i=1:5
    [x,fval,exitflag,output] = fminsearch( @ftsai_min, x, options, params_const, p3d, p2d)
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
    
    
    if size(angles,1)~=3
        angles=angles';
    end
    
    rotm = eul2rotm(angles);
    camParaCalib_new.R = rotm;

    camParaCalib_new.T = params(4:6)';
    
    
    
    T = [[camParaCalib_new.R params(4:6)']; [0 0 0 1]];
    Tinv = T^(-1);
    camParaCalib.Rinv = Tinv(1:3, 1:3);
    camParaCalib.Tinv = Tinv(1:3, 4);



    camParaCalib_new.f_eff = params(7);
    camParaCalib_new.k1 = params(8);
    camParaCalib_new.Noffw = params(9);
    camParaCalib_new.Noffh = params(10);
    camParaCalib_new.Rinv = Tinv(1:3, 1:3);
    camParaCalib_new.Tinv = Tinv(1:3, 4);
    

    
        Xproj = calibProj_Tsai(camParaCalib_new, p3d);
        err_img = Xproj - p2d;
        camParaCalib_new.err_x = sqrt(mean((err_img(:,1)).^2));
        camParaCalib_new.err_y = sqrt(mean((err_img(:,2)).^2));
        camParaCalib_new.err_t = sqrt((camParaCalib_new.err_x)^2 + (camParaCalib_new.err_y)^2);  
end