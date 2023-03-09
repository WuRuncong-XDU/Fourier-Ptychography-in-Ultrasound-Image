% -------------------------------------------------------------------------
% Paper:
% An Improved Fourier Ptychography Algorithm for Ultrasonic Array Imaging.
% Update 2023/03/09, Wu Runcong.
% -------------------------------------------------------------------------
clc, clear, close all;

%% Clear worksapce
AddFullPath;
if ~isempty(dir('file\PbestRecord.txt'))
    delete('file\PbestRecord.txt');
end
if ~isempty(dir('file\PgRecord.txt'))
    delete('file\PgRecord.txt');
end

%% Open parallel threads pool
fprintf("Parallel runing mode.\n");
if isempty(gcp('nocreate'))
    parpool('local');
else
    fprintf('CPU parallel pool has been opened.\n');
end

%% Set initial parameters and select optimate algorithm type.
% global ProbeType;
% ProbeType = 'HFL40';
OptimationType = 1;     % 1(PSO), 2(CSA), 3(MPA), 4(EOA), 5(GOA)
matFilename = {'HFL40_Simulation1.mat', 'HFL40_PigEye.mat', 'HFL40_Tungsten.mat', 'HFL40_Simulation2.mat'};	% Optional source data.
PreProcData(matFilename);

%% Iteration
% Target, unit [mm]
Des.resH = 0.010;
Des.resR = 0.010;
Des.contrast = -50;

% Population setting.
Dim = 3;
N = 25;
MaxIT = 50;
RunNo = 15;     % Only for EOA

switch OptimationType
    case 1
        % PSO-LDIW
        [BestFitness, BestPosition, PSOConvCurve] = PSO(N, MaxIT, Dim, Des);
    case 2
        % CSA
        [BestFitness, BestPosition, CSAConvCurve] = CSA(N, MaxIT, Dim, Des);
    case 3  
        % MPA
        [BestFitness, BestPosition, MPAConvCurve] = MPA(N, MaxIT, Dim, Des);
    case 4
        % EOA
        [BestFitness, BestPosition, EOAConvCurve] = EOA(N, MaxIT, Dim, RunNo, Des);
    case 5
        % GOA
        [BestFitness, BestPosition, GOAConvCurve] = GOA(N, MaxIT, Dim, Des);
    otherwise
        fprintf('Select a wrong algorithm type./n');
end

%% Use parameters after optimating
% [OutputParam, OutputImage, RefImage] = Ptychography30MHz(BestPosition);
% [OutputParam, OutputImage, RefImage] = Ptychography30MHzSimulation(BestPosition);
[OutputParam, OutputImage, RefImage] = Ptychography30MHzTungsten(BestPosition);

%% Save
% save('Output.mat', 'OutputImage', 'OutputParam', 'BestPosition', 'BestFitness');
