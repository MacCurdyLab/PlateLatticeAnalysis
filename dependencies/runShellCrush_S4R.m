function [simData] = runShellCrush_S4R(simStruct)

tic

hLatt = simStruct.h*(simStruct.nc(3)-1);
simStruct.vCrush = 100;           %[mm/s]

simStruct.dCrush = hLatt*simStruct.estrain;          %[mm]
simStruct.t_sim = simStruct.dCrush/simStruct.vCrush;

%build lattice geometry
[F,V]=makePlanarPlateLattice(simStruct.nc,simStruct.w,simStruct.h,simStruct.e);
V = V-mean(V);
A = sum(patchArea(F,V));
simData.A = A;

%compute lattice volume
wt = simStruct.wt;
LatVol = wt*A;
L = axisLim(V);
simData.L = L;
dims = diff(reshape(L,2,3));
simData.dims = dims;
TotVol = prod(dims);
simData.rho = LatVol/TotVol;

[F,V] = subQuad(F,V,simStruct.nSub);

noiz = simStruct.defectMag*wt;

V(:,1:2) = V(:,1:2)+noiz*(rand(size(V(:,1:2)))-0.5);

simData.F = F;
simData.V = V;

if simStruct.meshOnly
    return
end


fprintf('Mesh Complete with %1.2fk dof, rho = %f\n',numel(V)*1e-3,simData.rho)

fprintf('\nWriting INP File... \n')


%% INP WRITER

name = simStruct.name;
mkdir(name); 
cd(name);

fileID = fopen([name '.inp'],'w');

%% Write the Header
fprintf(fileID,'*Heading\n*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');

plateDim = [L(2)-L(1) L(4)-L(3)]*1.5;
[Fplate,Vplate]=quadPlate(plateDim,[5 5]);

boxThick = plateDim(1)/10;
[meshStruct]=hexMeshBox([plateDim boxThick],[5 5 1]);
Ebox=meshStruct.E;
Vbox=meshStruct.V;
Fbox=meshStruct.F;
Fbbox=meshStruct.Fb;

% cFigure; axisGeom;
% gpatch(F,V,'bw','k',1.0); hold on
% % gpatch(Fplate,Vplate+[0 0 L(5)],'gw','k',0.5);
% % gpatch(Fplate,Vplate+[0 0 L(6)],'gw','k',0.5);
% gpatch(Fbbox,Vbox+[0 0 L(6)+boxThick/2],'rw','k',0.5)
% gpatch(Fbbox,Vbox+[0 0 L(5)-boxThick/2],'rw','k',0.5)

%% Write the top plane
fprintf(fileID,'*Part, name=TOPPLATE\n*Node\n');
%node list
formatSpec = '%i,\t%f,\t%f,\t%f\n';
for i = 1:size(Vbox,1)
   line = [i Vbox(i,:)+[0 0 L(6)+boxThick/2]];
   fprintf(fileID,formatSpec,line);
end

%element list
fprintf(fileID,'*Element, type=C3D8R\n');
formatSpec = '%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i\n';
for i = 1:size(Ebox,1)
   line = [i Ebox(i,:)];
   fprintf(fileID,formatSpec,line);
end

%element set
fprintf(fileID,'*Elset, elset=allEl, generate\n');
fprintf(fileID,'1, \t%i,\n',size(Fplate,1));
fprintf(fileID,'*Solid Section, elset=allEl, material=SS_TPU\n');
fprintf(fileID,',\n');
fprintf(fileID,'*End Part\n');
% fprintf(fileID,'*Shell Section, elset=allEl, material=SS_TPU\n%1.5f, 5\n*End Part\n',0.1);

%% Write the bottom plane
fprintf(fileID,'*Part, name=BOTTOMPLATE\n*Node\n');
%node list
formatSpec = '%i,\t%f,\t%f,\t%f\n';
for i = 1:size(Vbox,1)
   line = [i Vbox(i,:)+[0 0 L(5)-boxThick/2]];
   fprintf(fileID,formatSpec,line);
end

%element list
fprintf(fileID,'*Element, type=C3D8R\n');
formatSpec = '%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i,\t%i\n';
for i = 1:size(Ebox,1)
   line = [i Ebox(i,:)];
   fprintf(fileID,formatSpec,line);
end

%element set
fprintf(fileID,'*Elset, elset=allEl, generate\n');
fprintf(fileID,'1, \t%i,\n',size(Fplate,1));
fprintf(fileID,'*Solid Section, elset=allEl, material=SS_TPU\n');
fprintf(fileID,',\n');
fprintf(fileID,'*End Part\n');
% fprintf(fileID,'*Shell Section, elset=allEl, material=SS_TPU\n%1.5f, 5\n*End Part\n',0.1);

%% Write the plate lattice
fprintf(fileID,'*Part, name=SQUAREBUCKLE\n*Node\n');

%node list
formatSpec = '%i,\t%f,\t%f,\t%f\n';
for i = 1:size(V,1)
   line = [i V(i,:)];
   fprintf(fileID,formatSpec,line);
end

%element list
fprintf(fileID,'*Element, type=S4R\n');
formatSpec = '%i,\t%i,\t%i,\t%i,\t%i\n';
for i = 1:size(F,1)
   line = [i F(i,:)];
   fprintf(fileID,formatSpec,line);
end

%element set
fprintf(fileID,'*Elset, elset=allEl, generate\n');
fprintf(fileID,'1, \t%i,\n',size(F,1));

fprintf(fileID,'*Shell Section, elset=allEl, material=SS_TPU\n%1.1f, 5\n*End Part\n',wt);

%% Assembly

fprintf(fileID,'*Assembly, name=Assembly  \n');
fprintf(fileID,'*Instance, name=SQUAREBUCKLE-1, part=SQUAREBUCKLE\n');
fprintf(fileID,'*End Instance  \n');
fprintf(fileID,'*Instance, name=TOPPLATE-1, part=TOPPLATE\n');
fprintf(fileID,'*End Instance  \n');
fprintf(fileID,'*Instance, name=BOTTOMPLATE-1, part=BOTTOMPLATE\n');
fprintf(fileID,'*End Instance  \n');
fprintf(fileID,'*Node\n');
fprintf(fileID,'      1,           0.,           0.,           %2.2f\n',L(5)*2);
fprintf(fileID,'*Nset, nset=m_Set-1\n');
fprintf(fileID,' 1,\n');
fprintf(fileID,' *Node\n');
fprintf(fileID,'      2,           0.,           0.,           %2.2f\n',L(6)*2);
fprintf(fileID,'*Nset, nset=m_Set-2\n');
fprintf(fileID,' 2,\n');

%couplings
fprintf(fileID,'*Elset, elset=_s_Surf-1_S5, internal, instance=TOPPLATE-1, generate\n');
fprintf(fileID,'  1,  25,   1\n');
fprintf(fileID,'*Surface, type=ELEMENT, name=s_Surf-1\n');
fprintf(fileID,'_s_Surf-1_S5, S5\n');
% fprintf(fileID,'*Elset, elset=_s_Surf-1_SNEG, internal, instance=TOPPLATE-1, generate\n');
% fprintf(fileID,'  1,  25,   1\n');
% fprintf(fileID,'*Surface, type=ELEMENT, name=s_Surf-1\n');
% fprintf(fileID,'_s_Surf-1_SNEG, SNEG\n');


fprintf(fileID,'*Elset, elset=_s_Surf-1_S6, internal, instance=BOTTOMPLATE-1, generate\n');
fprintf(fileID,'  1,  25,   1\n');
fprintf(fileID,'*Surface, type=ELEMENT, name=s_Surf-2\n');
fprintf(fileID,'_s_Surf-1_S6, S6\n');

% fprintf(fileID,'*Elset, elset=_s_Surf-2_SPOS, internal, instance=BOTTOMPLATE-1, generate\n');
% fprintf(fileID,'  1,  25,   1\n');
% fprintf(fileID,'*Surface, type=ELEMENT, name=s_Surf-2\n');
% fprintf(fileID,'_s_Surf-2_SPOS, SPOS\n');

fprintf(fileID,'*Coupling, constraint name=C1, ref node=m_Set-1, surface=s_Surf-1\n');
fprintf(fileID,'*Kinematic\n');

fprintf(fileID,'*Coupling, constraint name=C2, ref node=m_Set-2, surface=s_Surf-2\n');
fprintf(fileID,'*Kinematic\n');

fprintf(fileID,'*Nset, nset=allNodes, instance=SQUAREBUCKLE-1, generate\n');
fprintf(fileID,'\t1,\t%i,\t1\n',size(V,1));

fprintf(fileID,'*End Assembly\n');

fprintf(fileID,'*Material, name=SS_TPU\n*Density\n%2.2e,\n',simStruct.densityScale*1.21e-09);
fprintf(fileID,'*Hyperelastic, n=%i, ogden\n',length(simStruct.OgdenParams)/3);
fprintf(fileID,'%1.4f, ',simStruct.OgdenParams(1:end-1));
fprintf(fileID,'%1.4f\n',simStruct.OgdenParams(end));

if ~isempty(simStruct.g)
fprintf(fileID,'*Viscoelastic, time=PRONY\n');
for i = 1:size(simStruct.PronyParams,1)
fprintf(fileID,'  %1.4f,    %1.4f,  %1.4f\n',simStruct.PronyParams(i,:));
end
end

fprintf(fileID,'*Surface Interaction, name=conprops\n*Friction\n');
fprintf(fileID,'%2.2f,\n',simStruct.frictionCoeff);

fprintf(fileID,'*Surface Behavior, pressure-overclosure=HARD\n');
fprintf(fileID,'*Contact, op=NEW\n');
fprintf(fileID,'*Contact Inclusions, ALL ELEMENT BASED\n');
fprintf(fileID,'*Contact Property Assignment\n');
fprintf(fileID,' ,  , conprops\n');

fprintf(fileID,'*Step, name=Step-1, nlgeom=YES\n');
fprintf(fileID,'*Dynamic, Explicit\n');
fprintf(fileID,', %1.6f\n',simStruct.t_sim);

dt_target = simStruct.dt_target;
if dt_target ~= 0
fprintf(fileID,'*Variable Mass Scaling, dt=%2.6f, type=below min, frequency=5\n',dt_target);
end

fprintf(fileID,'*Bulk Viscosity\n');
fprintf(fileID,'0.12, 2.4\n');

fprintf(fileID,'*Boundary, type=VELOCITY\n');
fprintf(fileID,'m_Set-1, 1, 1\n');
fprintf(fileID,'m_Set-1, 2, 2\n'); 
fprintf(fileID,'m_Set-1, 3, 3, %2.6f\n',-simStruct.vCrush);
fprintf(fileID,'m_Set-1, 4, 4\n');
fprintf(fileID,'m_Set-1, 5, 5\n');
fprintf(fileID,'m_Set-1, 6, 6\n');

fprintf(fileID,'*Boundary\n');
fprintf(fileID,'m_Set-2, ENCASTRE\n');

fprintf(fileID,'*Restart, write, number interval=1, time marks=NO\n');
fprintf(fileID,'*FILE OUTPUT,NUMBER INTERVAL=');
fprintf(fileID,'%i\n',simStruct.nWriteFil);

fprintf(fileID,'*NODE FILE, NSET=m_Set-1\nRF, U\n');

fprintf(fileID,'*Output, field, variable=PRESELECT, number interval=%i\n',simStruct.nWriteODB);


fprintf(fileID,'*End Step');
fclose(fileID);

fprintf('Input File Complete... Dispatching Job\n')


%% Run it.
cmd_str = ['abaqus job=', name, ' input=', [name '.inp'] ' cpus=' num2str(simStruct.nCores) ' interactive'];
system(cmd_str);

% set up abaqus2Matlab
dir_path = pwd;
run('Documentation.m');
cd(dir_path);

cmd_str = ['abaqus ascfil job=', name];
system(cmd_str);

%% Extract the Outputs

try
Rec = Fil2str([name '.fin']);
catch
end

%extract time vector
try
out = Rec2000(Rec); 
simData.T = cell2mat(out(:,1)); 

out = Rec101(Rec);
U = out(:,4);
simData.U = U;

RF = stack3D(Rec104(Rec));
simData.RF = squeeze(RF(:,3,:))';

catch
    simData.T = []; 
    simData.RF = []; 
    simData.U = [];
end

cd ..
if simStruct.deleteMe
cmd_rmdir(name);
pause(1)
end

simData.runtime = toc;
simData.simStruct = simStruct;

save([simStruct.name '_results.mat'],"simData")
end
