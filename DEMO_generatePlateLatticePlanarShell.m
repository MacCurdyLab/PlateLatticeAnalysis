%% Demo: Shell Plate Lattice Generation
%Lawrence Smith | lasm4254@colorado.edu

simStruct.e = [0.75];      %[] geometric prebuckling (1 = perfectly orthogonal)
 simStruct.nc = [4 4 6]; %[] number of unit cells [x y z]
simStruct.h = 8;      %[mm] cell height      
simStruct.w = 12;       %[mm] cell edge length

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

%visualize plate lattice
cFigure; axisGeom; hold on
gpatch(F,V,'w','k',0,2) %bold outline
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
gdrawnow
