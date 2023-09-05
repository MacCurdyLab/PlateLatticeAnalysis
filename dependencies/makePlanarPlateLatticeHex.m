function [E, VE] = makePlanarPlateLatticeHex(nc,w,h,e,wt)

% cFigure; axisGeom;

n_x = nc(1); %number of unit cells, X
n_y = nc(2); %number of unit cells, Y
n_z = nc(3); %number of unit cells, Z

%determine the eccentricity vector
if length(e) == 1
    ee = e*ones(n_z,1);
elseif length(e) == 2
    ee = linspace(e(1),e(2),n_z);
else
    ee = interp1(linspace(1,n_z,length(e)),e,1:n_z);
end

V = []; F = []; E =[]; VE = [];

FF =[-4 -3 1 0; -3 -2 2 1; -2 -1 3 2; -1 -4 0 3];

flip1 = false;
flip2 = false;
flip3 = true;

for nx = 1:n_x

flip2 = ~flip2; 
flip3 = ~flip3;

if flip3
    flip1 = ~flip1;
end

for ny = 1:n_y

for nz = 1:n_z

    e = ee(nz);
    w2 = w/(1+e);
    m1 = (2*w2-w)/2;

    S1a = [[w/2 m1 -w/2 -m1]' [-m1 w/2 m1 -w/2]'];
    S1b = [[w/2 -m1 -w/2 m1]' [m1 w/2 -m1 -w/2]'];
    S2a = [[w/2-m1 0 -w/2+m1 0]' [0 w/2+m1 0 -w/2-m1]'];
    S2b = [[w/2+m1 0 -w/2-m1 0]' [0 w/2-m1 0 -w/2+m1]'];

    flip1 = ~flip1;
    Z = ones(4,1)*(nz-1)*h;

if flip2
    if flip1
    x = S2a(:,1); y = S2a(:,2);
    else
    x = S2b(:,1); y = S2b(:,2);
    end
else
    if flip1
    x = S1a(:,1); y = S1a(:,2);
    else
    x = S1b(:,1); y = S1b(:,2);
    end
end

X = x+w*(nx-1)/2;
Y = y+w*(ny-1);

if flip2
Y = Y+w/2;
end

%Let's not add a column if we are in the bad zone
if ~(iseven(nx) && ny==1)
    V = [V;[X Y Z]];
    if nz>1
        Vs = size(V,1)-3; 
        F = [F; Vs+FF;]; 
    end
end


end

if ~(iseven(nx) && ny==1)
    [E1,VE1]=quadThick(F,V,-1,wt/2,1);
    if ~isempty(E)
        [E,VE] = joinElementSets({E,E1},{VE,VE1});
    else
        E = E1; VE=VE1;
    end
end

% cla
% cFigure; axisGeom;
% gpatch(F,V,'bw',2,0.5); axis equal
% hold on
% gdrawnow

F = [];
V = [];

end

end

[E, VE] = mergeVertices(E,VE);

%findMinMaxVerticesXY

L = axisLim(VE);

Vout = find(VE(:,1) == L(1));
Vout = [Vout; find(VE(:,1) == L(2));];
Vout = [Vout; find(VE(:,2) == L(3));];
Vout = [Vout; find(VE(:,2) == L(4));];

Vout = unique(Vout);

FE=element2patch(E);
[indFree]=freeBoundaryPatch(FE);

AA = patchArea(FE,VE);

%we need to find the faces in the mesh that are on the exterior boundary.
%Here's the strategy. Start with a face that we know is on the boundary.
%Run face connectivity test to find all the adjacent faces. If any of these
%faces are not already on the list of keepers AND have an area that is
%greater than the average face area, they are added to the keepers list.
Keepers = find(any(Vout(1)==FE,2),1);
Seeded = Keepers;
Benched = [];

[C]=patchConnectivity(FE,VE);

FMap = C.face.face;

stopLoop = 0;
% cFigure; axisGeom
while ~stopLoop
%     cla
%     gpatch(FE(Keepers,:),VE,'bw') %plot all faces in set
%     gpatch(FE(Seeded,:),VE,'gw') %plot seed face
    new = FMap(Seeded,:);                   %all possible connected
    new = new(new>0 & ismember(new,indFree));   %cull faces not on boundary
    new = new(AA(new)>mean(AA));                %cull small faces
    new = new(~ismember(new,Keepers));          %cull faces already belonging to Keepers
%     gpatch(FE(new,:),VE,'rw')
%     gdrawnow();

    Keepers = [Keepers new];
    
    nBenched = new;
    
    if ~isempty(nBenched)
    nBenched(1) = [];
    Benched = [Benched nBenched];
    end
    
    if isempty(new)
        if isempty(Benched)
        stopLoop = 1;
        else
            Seeded = Benched(1);
            Benched(1) = [];
        end
    else
        Seeded = new(1);
    end

end

[E1,VE1]=quadThick(FE(Keepers,:),VE,-1,wt/2,1);
[E,VE] = joinElementSets({E,E1},{VE,VE1});
[E, VE] = mergeVertices(E,VE);

FE=element2patch(E);

[E,L_fixed]=fixNormalsOutward(E,VE,'s');

% cFigure; axisGeom
% gpatch(FE,VE,'gw','k',0.85)
% hold on
% plotV(VE(Vout,:),'r.','markersize',30)
% gdrawnow

end
