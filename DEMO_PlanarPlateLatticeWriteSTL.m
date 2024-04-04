%% DEMO - Create an STL of a Plate Lattice

clear; clc; close all

plotme = false; %should we generate a plot of the lattice?

nc = [7 4 7];   %number of cells in x, y, and z
%Note - we have not tested the robustness to different combinations of cell
%numbers in x,y,z. You might need to keep the number of cells in X odd and
%the number of cells in Y even for example.

w = 12;         %[mm] cell width
h = 8;          %[mm] cell height
e = 0.75;       %[] eccentricity (1-prebuckling factor)
wt = 1;         %[mm] wall thickness

%create a filename for the STL
if length(e) == 1
name = sprintf('w%i_h%i_e%i_wt0%i',[w h e*100 wt*100]);
else
name = sprintf('w%i_h%i_e%i-%i_wt0%i',[w h e(1)*100 e(2)*100 wt*100]);
end

% sprinkle some fairy dust on the wall thickness
magic_scale = 0.8485/0.6;
wt = wt*magic_scale;

%Compute
[E, VE] = makePlanarPlateLatticeHex(nc,w,h,e,wt);

%Find the free boundary faces
[F]=element2patch(E,'hex8');
if plotme
cFigure; axisGeom
gpatch(F,VE,'gw','k',0.85); hold on
end
[indFree]=freeBoundaryPatch(F);

%Find which patches are large enough
AA = patchArea(F,VE);
largeEnough = AA'>mean(AA);
isExterior = ismember(1:size(F,1),indFree);

%now lets perform a calculation of the actual wall thickness. We need to
%find two faces located on the exterior of the lattice and compute the
%distance between their centers.
C = patchCentre(F,VE);
C(~largeEnough) = -inf+C(~largeEnough);
C(~isExterior) = -inf+C(~isExterior);
[~,i] = sort(sum(C,2),'descend');

a = C(i(1),:);
b = C(i(2),:);
if plotme
plotV(a,'r.','markersize',20)
plotV(b,'b.','markersize',20)
end
dist = sqrt(sum((a-b).^2));

fprintf('\nActual wall thickness is %2.4f mm\n',dist)

[F,V] = patchCleanUnused(F(indFree,:),VE);

%subdivide the plate lattice faces
[F,V] = subQuad(F,V,4);
[Ft,Vt]=quad2tri(F,V,'f');

%reverse the normals
Ft = Ft(:,[1 3 2]);

%write an STL for the lattice
stlwrite(triangulation(Ft,Vt),['SOLID_' name '.stl'])


