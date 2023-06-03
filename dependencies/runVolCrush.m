function simData = runVolCrush(simStruct)

% I'm envisioning something that is more generic and can take in a file for
% the NC and CM matricies separately instead of being tied to a single
% mesh.

tic
name = simStruct.name;
mkdir(name)
cd(name)


CM = importdata(['dep/' simStruct.CMFile]);
[F,~]=element2faces(CM(:,2:end),[]);
NC = importdata(['dep/' simStruct.NCFile]);
simData.F = F;
simData.E = CM(:,2:end);
simData.V = NC(:,2:end);

if simStruct.meshOnly
    cd ..
    return
end

fileID = fopen([name '.inp'],'w');

fprintf(fileID,'*Heading\n');
fprintf(fileID,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');

fid = fopen('../dep/addPlane.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);

fprintf(fileID,'*Part, name=PlateLattice\n');
fprintf(fileID,'*Node\n');

fprintf(fileID,'\t%d,\t %1.3f,\t %1.3f,\t %1.3f,\t\n',NC');

fprintf(fileID,'*Element, type=C3D8R\n');

fprintf(fileID,'\t%d,\t %d,\t %d,\t %d,\t %d,\t %d,\t %d,\t %d,\t %d,\t\n',CM');

fprintf(fileID,'*Nset, nset=PL_section, generate\n');
fprintf(fileID,'      1,  %d,       1\n',size(NC,1));
fprintf(fileID,'*Elset, elset=PL_section, generate\n');
fprintf(fileID,'      1,  %d,       1\n',size(CM,1));
fprintf(fileID,'*Solid Section, elset=PL_section, material=SS_TPU\n');
fprintf(fileID,',\n');
fprintf(fileID,'*End Part\n');

fid = fopen('../dep/AssemblySets.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);

fprintf(fileID,'*End Assembly\n**\n');
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

fprintf(fileID,'*Surface Interaction, name=ContProps\n');
fprintf(fileID,'*Friction\n');
fprintf(fileID,'%2.2f,\n',simStruct.frictionCoeff);
fprintf(fileID,'*Surface Behavior, pressure-overclosure=HARD\n');

fprintf(fileID,'*Contact, op=NEW\n');
fprintf(fileID,'*Contact Inclusions, ALL EXTERIOR\n');
fprintf(fileID,'*Contact Property Assignment\n');
fprintf(fileID,' ,  , ContProps\n');

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
fprintf(fileID,'m_Set-3, ENCASTRE\n');

fprintf(fileID,'*Restart, write, number interval=1, time marks=NO\n');
fprintf(fileID,'*FILE OUTPUT,NUMBER INTERVAL=');
fprintf(fileID,'%i\n',simStruct.nWriteFil);

fprintf(fileID,'*NODE FILE, NSET=m_Set-1\nRF, U\n');

fprintf(fileID,'*Output, field, variable=PRESELECT, number interval=%i\n',simStruct.nWriteODB);

fprintf(fileID,'*End Step');
fclose(fileID);

fprintf('Input File Complete... Dispatching Job\n')


cmd_str = ['abaqus job=', name, ' input=', [name '.inp'] ' cpus=' num2str(simStruct.nCores) ' interactive'];
system(cmd_str);
fprintf('\nSimulation Completed in %1.1f s \n',toc)

% set up abaqus2Matlab
dir_path = pwd;
run('Documentation.m');
cd(dir_path);

cmd_str = ['abaqus ascfil job=', name];
system(cmd_str);
fprintf('\nFile Conversion Completed in %1.1f s \n',toc)

try
Rec = Fil2str([name '.fin']);
catch
end

try
%extract time vector
out = Rec2000(Rec); 
T = cell2mat(out(:,1));

out = Rec104(Rec); 
RF = out(:,4);

out = Rec101(Rec);
U = out(:,4);

simData.T = T;
simData.RF = RF;
simData.U = U;

catch
    simData.T = []; 
    simData.RF = []; 
    simData.U = [];
end

simData.runtime = toc;
simData.simStruct = simStruct;

cd ..
if simStruct.deleteMe
cmd_rmdir(name);
pause(1)
end


end