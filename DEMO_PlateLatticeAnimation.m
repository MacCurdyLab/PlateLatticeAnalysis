%% Demo: Shell Plate Lattice Animation Generation
%Lawrence Smith | lasm4254@colorado.edu

clear; clc; close all;

simStruct.e = [0.75];      %[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.nc = [8 8 8]; %[] number of unit cells [x y z]
simStruct.h = 8;      %[mm] cell height      
simStruct.w = 12;       %[mm] cell edge length

%visualize plate lattice
cFigure; axisGeom; hold on
set(gca,'visible','off')

x = linspace(0,1,10);

k = 0;

%% Smoothly Vary one Way
for i = linspace(0,1,25)
k = k+1;

%[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.e = 1-0.5*i;      

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

cla
gpatch(F,V,'w','k',0,2) %bold outline
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
campos([-171.5866 -352.7406  111.1608])
gdrawnow

img = getframe();
img.cdata = imresize(img.cdata,[800,1200]);
FF(k) = img;

end
%% Smoothly Vary The other Way
for i = linspace(0,1,25)
k = k+1;

%[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.e = 0.5+0.5*i;      

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

cla
gpatch(F,V,'w','k',0,2) %bold outline
axis manual
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
campos([-171.5866 -352.7406  111.1608])
gdrawnow

img = getframe();
img.cdata = imresize(img.cdata,[800,1200]);
FF(k) = img;

end

%% Introduce a little buckle at the bottom
for i = linspace(0,1,20)
k = k+1;

%[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.e = 1-i*0.5*rescale(normpdf(x,1,0.15));      

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

cla
gpatch(F,V,'w','k',0,2) %bold outline
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
campos([-171.5866 -352.7406  111.1608])
gdrawnow

img = getframe();
img.cdata = imresize(img.cdata,[800,1200]);
FF(k) = img;

end

%% Sweep the buckle upwards
for i = linspace(0,1,60)
k = k+1;

%[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.e = 1-0.5*rescale(normpdf(x,1-i,0.15));      

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

cla
gpatch(F,V,'w','k',0,2) %bold outline
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
campos([-171.5866 -352.7406  111.1608])
gdrawnow

img = getframe();
img.cdata = imresize(img.cdata,[800,1200]);
FF(k) = img;

end

%% Sweep the buckle upwards
for i = linspace(0,1,60)
k = k+1;

%[] geometric prebuckling (1 = perfectly orthogonal)
simStruct.e = 1-0.5*rescale(normpdf(x,i,0.15));      

%generate plate lattice
[F,V] = makePlanarPlateLatticeShell(simStruct);

cla
gpatch(F,V,'w','k',0,2) %bold outline
[F,V] = subQuad(F,V,2);
gpatch(F,V,'w','k',1,1) %face hatching
campos([-171.5866 -352.7406  111.1608])
gdrawnow

img = getframe();
img.cdata = imresize(img.cdata,[800,1200]);
FF(k) = img;

end


v = VideoWriter('NewVid.avi','Uncompressed AVI');
v.FrameRate = 20;
open(v)
writeVideo(v,FF)
close(v)