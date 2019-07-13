function [in,S] =particledetector(xcenter,Nx,ycenter,Ny,zcenter,Nz,particlespos)
    dx =[xcenter-Nx/2 xcenter+Nx/2]; 
    dy =[ycenter-Ny/2 ycenter+Ny/2]; 
    dz =[zcenter-Nz/2 zcenter+Nz/2]; 
    [x,y,z] = meshgrid(dx,dy,dz);  % A cube
    x = [x(:);xcenter];
    y = [y(:);ycenter];
    z = [z(:);zcenter];
    % [x,y,z] are corners of a cube plus the center.
    X = [x(:) y(:) z(:)];
    %Cubic 3D triangulation of each subvolume
    Tes = delaunayTriangulation(X); 
    %Define faces and vertices as inputs for inpolyhedron
    [S.faces, S.vertices]=freeBoundary(Tes); 
    % Logic list that particles exist in each subvolume
    in = inpolyhedron(S,particlespos); 

    