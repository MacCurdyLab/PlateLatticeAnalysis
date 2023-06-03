%% Demo: Run Impact Simulation 
% Lawrence Smith | lasm4254@colorado.edu

clear; clc; close all

%% Set up Simulation

%Geometry Settings
simStruct.e = [0.6 1.0];        %[] eccentricity (1 = perfectly orthogonal)
simStruct.nc = [4 4 6];         %[] number of unit cells [x y z]
simStruct.h = 8;                %[mm] cell height      
simStruct.w = 12;               %[mm] cell edge length
simStruct.wt = 0.6;             %[mm] PL wall thickness
simStruct.defectMag = 0.0;      %[] perturbation factor for buckling analysis, keep at 0.0
simStruct.nSub = 3;             %[] Number of subdivisions for each face, usually try 2 or 3 

%Impact Settings
simStruct.vImpact = 2500;       %[mm/s] impact velocity
simStruct.mImpact = 5.0;        %[kg]   impact mass
simStruct.frictionCoeff=0.75;   %[] global friction coefficien

%Simulation Settings
simStruct.name = 'testSim';
simStruct.nWriteFil = 200;      %[] this is the number of force-displacement probe results to write 
simStruct.nWriteODB = 40;       %[] the number of FULL-FIELD results (adds up fast)
simStruct.t_sim = 0.055;        %[s] total simulation time (usually 30-50ms)
simStruct.nCores = 10;          %[] number of cores to divide the simulation over  
simStruct.deleteMe = false;     %[] should we delete the simulation files?
simStruct.dt_target = 0;        %[] optional accelerator (5e-5 seems good, keep at 0 for impact)
simStruct.densityScale = 1;     %[] density override, keep at 1.0

%Material Properties: Ogden Hyperelasticity
mu =    [2.816490181E-02 7.81250824];
alpha = [4.24795953  -2.07903716];

%Material Properties: Prony Viscoelasticity
simStruct.tau = [0.001  0.0100   0.1000 ];
simStruct.g =   [0.5039   0.1863   0.0181]';    

%Store Material Properties in  simStruct
OgdenParams = [mu(:) alpha(:)]';
DOgden = zeros(size(mu));
simStruct.OgdenParams = [OgdenParams(:); DOgden(:)];
simStruct.PronyParams = [simStruct.g(:) 0*simStruct.g(:) simStruct.tau(:)];

%% Run Shell Sim
simDataShell = runShellImpact_S4R(simStruct);

%plot impact force vs. time
plot(simDataShell.T*1e3,simDataShell.RF','k-','linewidth',2,...
'displayname',sprintf('Shell FEA, t_{sim} = %1.1e s',simDataShell.runtime));
xlabel('Time [ms]')
ylabel('Impact Force [N]')
set(gca,'fontname','georgia','fontsize',14)
legend('location','southeast')
drawnow

