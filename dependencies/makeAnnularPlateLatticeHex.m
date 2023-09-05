function [E, VE] = makeAnnularPlateLatticeHex(nc,r,e,wt,varargin)

% cFigure; axisGeom;

n_c = nc(1); %number of unit cell pairs to tile around circumference
n_r = nc(2); %number of radial rings
n_h = nc(3); %number of vertical stacks
ro = r(1);
ri = r(2);

if isempty(varargin)
    plotme = false;
else
    plotme = varargin{1};
end



if plotme
cFigure; axisGeom
end

%determine the eccentricity vector
if length(e) == 1
    ee = e*ones(n_r+1,1);
elseif length(e) == 2
    ee = linspace(e(1),e(2),n_r+1);
else
    ee = interp1(linspace(1,n_r+1,length(e)),e,1:n_r+1);
end

%determine wall thickness
if length(wt)>1
    %we will perform a nudge so that the outer wt is equal to wt(2)
    performNudge = true;
    nudgeMag = (wt(2)-wt(1))/2;
    wt = wt(1);
else
    performNudge = false;
end

h = (ro-ri)/n_r;
w = pi*(ri+h*n_r/2)/n_c;

V = []; F = []; E = []; VE = [];

theta = 0;

FF =[-4 -3 1 0; -3 -2 2 1; -2 -1 3 2; -1 -4 0 3];
FF2 =[-4 -3 1 0; -3 -2 2 1; -2 -1 3 2; -1 -4 0 3]; 

flip1 = false;
flip2 = false;
flip3 = true;

for hx = 1:n_h*2

flip2 = ~flip2; 
flip3 = ~flip3;

if flip3
    flip1 = ~flip1;
end

for cx = 1:n_c*2

theta = (cx-1)*pi/n_c;

for rx = 1:n_r+1

    e = ee(rx);
    w2 = w/(1+e);
    m1 = (2*w2-w)/2;

    S1a = [[w/2 m1 -w/2 -m1]' [-m1 w/2 m1 -w/2]'];
    S1b = [[w/2 -m1 -w/2 m1]' [m1 w/2 -m1 -w/2]'];
    S2a = [[w/2-m1 0 -w/2+m1 0]' [0 w/2+m1 0 -w/2-m1]'];
    S2b = [[w/2+m1 0 -w/2-m1 0]' [0 w/2-m1 0 -w/2+m1]'];

    flip1 = ~flip1;
    R = ri + (rx-1)*h;
    fs = R/(0.5*(ri+ro));

if flip2
    if flip1
    x = S2a(:,1); y = S2a(:,2);
    else
    x = S2b(:,1); y = S2b(:,2);
    end
    x = x+w/2;
else
    if flip1
    x = S1a(:,1); y = S1a(:,2);
    else
    x = S1b(:,1); y = S1b(:,2);
    end
end

y = y+w/2*(hx-1);

[X,Y,Z] = pol2cart(theta+fs*x./R,R*ones(4,1),y);

V = [V;[X Y Z]];

if rx>1
    Vs = size(V,1)-3; 
    F = [F; Vs+FF;]; 
end

end

[E1,VE1]=quadThick(F,V,-1,wt/2,1);
% cFigure; axisGeom;
% gpatch(F,V,'bw','b',0.2); hold on
% plotV(VE1(1:size(V,1),:),'.','markersize',20,'color','r')
% plotV(VE1(size(V,1)+1:end,:),'.','markersize',20,'color','g')

% FQ = element2patch(E1,'hex8');
% gpatch(FQ,VE1,'gw','k',0.15); hold on
if performNudge
for i = 1:size(V,1)
    a = VE1(i,:);
    b = VE1(i+size(V,1),:);
    nudgeDist = sqrt(sum((b-a).^2,2));
    nudgeDir = (b-a)/nudgeDist;
    
    [~,mR]=cart2pol(b(1),b(2));
    nR = (mR-r(2))/(r(1)-r(2));

    %perform nudge
%     plotV(VE1(i+size(V,1),:),'.','markersize',20,'color','r')
    VE1(i+size(V,1),:) = b+nudgeDir*nudgeMag*nR;
%     plotV(VE1(i+size(V,1),:),'.','markersize',20,'color','b')
    
end
end
% FQ = element2patch(E1,'hex8');
% gpatch(FQ,VE1,'gw','k',0.15); hold on


if ~isempty(E)
    [E,VE] = joinElementSets({E,E1},{VE,VE1});
else
    E = E1; VE=VE1;
end

[F2]=element2patch(E,'hex8');
% cla
% gpatch(F2,VE,'gw','k',0.15); hold on

F = [];
V = [];

VEf = VE;
FEf = F2;

end

end

% close all

[E, VE] = mergeVertices(E,VE);

%findMinMaxVerticesXY
L = axisLim(VE);

%% Extrude a bit on the upper face
Vtop = find(VE(:,3) == L(6));

FE=element2patch(E);
[indFree]=freeBoundaryPatch(FE);
AA = patchArea(FE,VE);

%we need to find the faces in the mesh that are on the exterior boundary.
%Here's the strategy. Start with a face that we know is on the boundary.
%Run face connectivity test to find all the adjacent faces. If any of these
%faces are not already on the list of keepers AND have an area that is
%greater than the average face area, they are added to the keepers list.
Keepers = find(any(Vtop(1)==FE,2));
[~,ki] = sort(AA(Keepers),'descend');
Keepers = Keepers(ki(1));
Seeded = Keepers;
Benched = [];

[C]=patchConnectivity(FE,VE);

FMap = C.face.face;

stopLoop = 0;
% cFigure; axisGeom
while ~stopLoop
    if plotme
    cla
    gpatch(FEf,VEf,'gw','k',0.15); hold on
    gpatch(FE(Keepers,:),VE,'bw') %plot all faces in set
    gpatch(FE(Seeded,:),VE,'gw') %plot seed face
    end

    new = FMap(Seeded,:);                   %all possible connected
    new = new(new>0 & ismember(new,indFree));   %cull faces not on boundary
    new = new(AA(new)>mean(AA));                %cull small faces
    new = new(~ismember(new,Keepers));          %cull faces already belonging to Keepers
    
    if plotme
    gpatch(FE(new,:),VE,'rw')
    gdrawnow();
    end

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

[Fsmall,Vsmall] = patchCleanUnused(FE(Keepers,:),VE);

[E1,VE1]=quadThick(Fsmall,Vsmall,-1,wt/2,1);

if performNudge
for i = 1:size(Vsmall,1)
    a = VE1(i,:);
    b = VE1(i+size(Vsmall,1),:);
    nudgeDist = sqrt(sum((b-a).^2,2));
    nudgeDir = (b-a)/nudgeDist;
    
    [~,mR]=cart2pol(b(1),b(2));
    nR = (mR-r(2))/(r(1)-r(2));

    %perform nudge
    VE1(i+size(Vsmall,1),:) = b+nudgeDir*nudgeMag*nR;    
end
end

FlatBottom = true;

[FQ]=element2patch(E1,'hex8');

if plotme
cFigure; axisGeom
gpatch(FQ,VE1,'bw','k',0.55); hold on 
end

[~,si] = sort(VE1(size(Vsmall,1)+1:end,3),'descend');

si = si(1:size(Vsmall,1)/2);

if plotme
plotV(VE1(si,:),'.','markersize',20,'color','r')
end

if FlatBottom
VE1(si,3) = mean(VE1(si,3));
end

[E,VE] = joinElementSets({E,E1},{VE,VE1});
[E, VE] = mergeVertices(E,VE);

%% NOW REPEAT FOR LOWER FACE
Vtop = find(VE(:,3) == L(5));

FE=element2patch(E);
[indFree]=freeBoundaryPatch(FE);
AA = patchArea(FE,VE);

%we need to find the faces in the mesh that are on the exterior boundary.
%Here's the strategy. Start with a face that we know is on the boundary.
%Run face connectivity test to find all the adjacent faces. If any of these
%faces are not already on the list of keepers AND have an area that is
%greater than the average face area, they are added to the keepers list.
Keepers = find(any(Vtop(1)==FE,2));
[~,ki] = sort(AA(Keepers),'descend');
Keepers = Keepers(ki(1));
Seeded = Keepers;
Benched = [];

[C]=patchConnectivity(FE,VE);

FMap = C.face.face;

stopLoop = 0;
% cFigure; axisGeom
while ~stopLoop
    if plotme
    cla
    gpatch(FEf,VEf,'gw','k',0.15); hold on
    gpatch(FE(Keepers,:),VE,'bw') %plot all faces in set
    gpatch(FE(Seeded,:),VE,'gw') %plot seed face
    end
    new = FMap(Seeded,:);                   %all possible connected
    new = new(new>0 & ismember(new,indFree));   %cull faces not on boundary
    new = new(AA(new)>mean(AA));                %cull small faces
    new = new(~ismember(new,Keepers));          %cull faces already belonging to Keepers
    if plotme
    gpatch(FE(new,:),VE,'rw')
    gdrawnow();
    end

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

[Fsmall,Vsmall] = patchCleanUnused(FE(Keepers,:),VE);

[E1,VE1]=quadThick(Fsmall,Vsmall,-1,wt/2,1);

if performNudge
for i = 1:size(Vsmall,1)
    a = VE1(i,:);
    b = VE1(i+size(Vsmall,1),:);
    nudgeDist = sqrt(sum((b-a).^2,2));
    nudgeDir = (b-a)/nudgeDist;
    
    [~,mR]=cart2pol(b(1),b(2));
    nR = (mR-r(2))/(r(1)-r(2));

    %perform nudge
    VE1(i+size(Vsmall,1),:) = b+nudgeDir*nudgeMag*nR;    
end
end


[E,VE] = joinElementSets({E,E1},{VE,VE1});
[E, VE] = mergeVertices(E,VE);




end
