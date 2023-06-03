function [simData] = runVolImpact(simStruct)

tic

name = simStruct.name;
mkdir(name)
cd(name)

fileID = fopen([name '.inp'],'w');

fid = fopen('../dep/p1.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);

fprintf(fileID,'%1.5f,\n',simStruct.mImpact*1e-3); 

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

fid = fopen('../dep/p2.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);

fprintf(fileID,'-%1.2f',simStruct.vImpact); 

fid = fopen('../dep/p3.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);

fprintf(fileID,' %1.5f\n',simStruct.t_sim);

fid = fopen('../dep/p4.txt','r');
writeme = fread(fid, '*char'); 
fclose(fid);
fprintf(fileID,writeme);


fprintf(fileID,'%i\n',simStruct.nWriteODB);
fprintf(fileID,'*OUTPUT, HISTORY, TIME INTERVAL=1.0E-4\n');
fprintf(fileID,'*ENERGY OUTPUT\n');
fprintf(fileID,'ALLKE, ALLSE, ALLWK, ALLIE, ALLVD, ALLAE, ETOTAL\n');
fprintf(fileID,'*End Step\n');


fclose(fileID);

tic
cmd_str = ['abaqus job=', name, ' input=', [name '.inp'] ' cpus=' num2str(simStruct.nCores) ' interactive'];
system(cmd_str);
fprintf('\nSimulation Completed in %1.1f s \n',toc)

% set up abaqus2Matlab
dir_path = pwd;
run('Documentation.m');
cd(dir_path);

cmd_str = ['abaqus ascfil job=', name];
system(cmd_str);

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
catch
    T = []; RF = []; U = [];
end


cd ..
if simStruct.deleteMe
cmd_rmdir(name);
pause(1)
end

simData.T = T;
simData.RF = RF;
simData.U = U;
simData.runtime = toc;
simData.simStruct = simStruct;

end