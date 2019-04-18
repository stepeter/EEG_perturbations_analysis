%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_rotationsDenoise/Code/'));
    startUpEEGLAB('close');
end
close all;
%Specify parameters
subjnum=8; %subject number
RootData='/media/stepeter/Local_Data/VR_rotationsDenoise/Data/';

subjcode = ['WMISM_' num2str(subjnum)];
group='all'; %'mismatches'; %'pull';
subj_basepath = [RootData subjcode filesep]; % Top level subject folder
EEGsets_outpath = [subj_basepath 'EEG_sets']; % Where you want to save your .set files
ICA_path = [subj_basepath 'ICA_Stuff'];
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];

%Create directories (if necessary)
if ~exist(ICA_path)
    mkdir(ICA_path)
end
if ~exist(ICA_path_out)
    mkdir(ICA_path_out)
end
%% Load in .set file
cd(subj_basepath)
global ALLEEG EEG CURRENTSET ALLCOM;
EEG=pop_loadset('filename', ['Merge_' group '.set'], 'filepath', EEGsets_outpath);
numChansOrig=EEG.nbchan;

%Fix urchan issue with externals
if EEG.nbchan==133
    for i=129:133
        EEG.chanlocs(1,i).urchan=i;
    end
end

EEG_orig=EEG;

cd(ICA_path_out) %cd to save rejection log files in ICA_path_out
if exist([ICA_path_out filesep group '_ICA_PREP_1.mat']) == 2
    krej = input('You have already rejected channels once, keep this rejection? (yes = 0, no = 1): ');
    
    if krej == 0
        load([ICA_path_out filesep group '_ICA_PREP_1.mat']);
    elseif krej == 1
        
        % Reject bad channels (head and externals separately)
        if EEG.nbchan==128
            [ICA_STRUCT] = chan_rej_hjh(EEG); %With all scalp electrodes
            % Add good externals to ICA_STRUCT.good_chans
            ICA_STRUCT.good_cap = ICA_STRUCT.good_chans; % Store non-externals into 'good_cap'
            ICA_STRUCT.good_ext = []; % No externals
        elseif EEG.nbchan==133
            [ICA_STRUCT] = chan_rej_hjh(EEG,129:133);
            [ICA_STRUCT_ext] = chan_rej_hjh(EEG,1:128); %just with externals
            ICA_STRUCT.good_cap = ICA_STRUCT.good_chans; % Store non-externals into 'good_cap'
            ICA_STRUCT.good_ext=ICA_STRUCT_ext.good_chans;
        end
        
            
    %updating ICA structure
    ICA_STRUCT.filename = [group '_ICA_PREP_1.mat'];
    ICA_STRUCT.filepath = ICA_path_out;
    end
    
else
    display('No previous channel rejection found...')
    % Reject bad channels (head and externals separately)
    if EEG.nbchan==128
        [ICA_STRUCT] = chan_rej_hjh(EEG); %With all scalp electrodes
        % Add good externals to ICA_STRUCT.good_chans
        ICA_STRUCT.good_cap = ICA_STRUCT.good_chans; % Store non-externals into 'good_cap'
        ICA_STRUCT.good_ext = []; % No externals
    elseif EEG.nbchan==133
        [ICA_STRUCT] = chan_rej_hjh(EEG,129:133);
        [ICA_STRUCT_ext] = chan_rej_hjh(EEG,1:128); %just with externals
        ICA_STRUCT.good_cap = ICA_STRUCT.good_chans; % Store non-externals into 'good_cap'
        ICA_STRUCT.good_ext=ICA_STRUCT_ext.good_chans;
    end       
            
    %updating ICA structure
    ICA_STRUCT.filename = [group '_ICA_PREP_1.mat'];
    ICA_STRUCT.filepath = ICA_path_out;
end
cd(subj_basepath)
%% Interpolate and re-reference
% Interpolate bad channels
EEG = update_EEG(EEG, ICA_STRUCT);
EEG=pop_interp(EEG,EEG_orig.chanlocs,'spherical');

% Save ICA_STRUCT with just bad channels rejected
save([ICA_path_out filesep group '_ICA_PREP_1'], 'ICA_STRUCT')

% %Re-reference externals if needed (do this before CAR)
% if EEG.nbchan==133
%     neckAve=(EEG.data(129,:)+EEG.data(130,:))/2;
%     EEG.data(129,:)=EEG.data(129,:)-neckAve;
%     EEG.data(130,:)=EEG.data(130,:)-neckAve;
%     eyeAve=(EEG.data(131,:)+EEG.data(132,:)+EEG.data(133,:))/3;
%     EEG.data(131,:)=EEG.data(131,:)-eyeAve;
%     EEG.data(132,:)=EEG.data(132,:)-eyeAve;
%     EEG.data(133,:)=EEG.data(133,:)-eyeAve;
% end

%Reference to a common average and plot spectra
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:) = zeros(1, EEG.pnts);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select( EEG,'nochannel',{'initialReference'});
Vfig=figure;
pop_spectopo(EEG, 1, [0 EEG.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
saveas(Vfig,[EEGsets_outpath filesep 'PowerSpectrum' group '_channel_rej'],'fig');

fprintf('Rejected %d channels\n',numChansOrig-length(ICA_STRUCT.good_chans))
disp('Channel Rejection Done :)');


%Save EEG set
pop_saveset(EEG, 'filename', ['Merge_' group '_CAR'], 'filepath', EEGsets_outpath,'savemode','twofiles');
%% Frame rejection

%Calculate number of frames needed to keep for k values of 60 and 30
target_nframes = length(ICA_STRUCT.good_chans)^2 * 60 / EEG.pnts; 
fprintf('You need to keep %.4f frames to achieve k=60\n',target_nframes)

target_nframes2 = length(ICA_STRUCT.good_chans)^2 * 30 / EEG.pnts; 
fprintf('You need to keep %.4f frames to achieve k=30\n',target_nframes2)

% Reject bad frames
[ICA_STRUCT] = frame_rej_hjh(EEG,ICA_STRUCT);

%% update ICA Struct
%Rereference
ICA_STRUCT.good_cap=1:EEG.nbchan; ICA_STRUCT.good_chans=1:EEG.nbchan;
ICA_STRUCT.ref = 'average'; 
EEG = update_EEG(EEG, ICA_STRUCT);

%Remove bad frames
EEG_temp=EEG;
EEG_temp.data(:,ICA_STRUCT.rej_frame_idx) = [];
EEG_temp.times = EEG_temp.times(1:size(EEG_temp.data,2));

%Plot power spectrum for frame rejected data
Vfig=figure;
pop_spectopo(EEG_temp, 1, [0 EEG_temp.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
saveas(Vfig,[EEGsets_outpath filesep 'PowerSpectrum' group '_frame_rej'],'fig');

%Plot final channel distribution after rejections
FigCh = figure;
topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
saveas(FigCh,[ICA_path_out filesep 'Channels_' subjcode '_' group],'fig') %Tsk(TaskType).name

%Save ICA_STRUCT with bad channels and bad frames
ICA_STRUCT.filename = [group '_ICA_PREP_2.mat'];
ICA_STRUCT.filepath = ICA_path_out;
save([ICA_path_out filesep group '_ICA_PREP_2'], 'ICA_STRUCT')


%% SAVE FLOAT FILE FOR ICA

% Save float file for running ICA
tmpdata = EEG.data;
tmpdata(:,ICA_STRUCT.rej_frame_idx) = []; %remove bad frames
disp('Converting data to double...');
tmpdata = double(tmpdata);
filename = [subjcode group '.fdt'];
data2write = tmpdata;
tmpdatafile = fullfile(ICA_path,filename);
disp(['Writing data file: ' tmpdatafile]);
floatwrite(data2write, tmpdatafile);
% %Write file to data3 as well
% tmpdatafile3 = fullfile(ICA_path3,filename);
% disp(['Writing data file: ' tmpdatafile3]);
% floatwrite(data2write, tmpdatafile3);
num_chans = size(tmpdata,1); % Number of channels used
disp(['Number of channels: ' num2str(num_chans)]);
num_frames = size(tmpdata,2); % Number of frames used
disp(['Number of frames: ' num2str(num_frames)]);
k = num_frames / num_chans^2; % K value achieved
disp(['Final k value: ' num2str(k)]);

%% Define parameter file for ICA

% RootCluster = '/sshfs/kin-hnl-frontend-stepeter/share/data4/stepeter/Projects/WMISM/Data';
RootFlux = '/ufrc/dferris/s.peterson/Projects/VR_rotations/Data';
% pGran = [RootCluster '/' subjcode '/ICA_Stuff/'];
pFlux = [RootFlux '/' subjcode '/ICA_Stuff/'];

% % Save CUDAICA parameter file to subject folder
% saveCUDAICAparamfile(ICA_path,pGran,filename,num_chans,num_frames,group);
% %Create CUDAICA folder
% cudaicaPath=[subj_basepath 'ICA_Stuff' filesep 'files_CUDAICA' group];
% if ~exist(cudaicaPath)
%     mkdir(cudaicaPath)
% end

% Save AMICA .param file to subject folder
maxthreads=4; %24
doPCA=1; numPCs=90;
if doPCA==1
    k = num_frames / numPCs^2; % K value achieved
    disp(['Final k value (with PCA): ' num2str(k)]);
end
saveAMICAparamfile(ICA_path,pFlux,filename,num_chans,num_frames,doPCA,numPCs,maxthreads,group);
saveAMICAparamfile(ICA_path3,pFlux,filename,num_chans,num_frames,doPCA,numPCs,maxthreads,group);

% Save qsub parameter file to subject folder (for running on Flux)
email='stepeter@umich.edu';
account='ferrisdp_flux';
procs=8; mem='4000mb'; walltime='12:00:00';
saveQSUBfile(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group);
saveQSUBfile(ICA_path3,pFlux,filename,email,account,procs,mem,walltime,group);

disp('Step 2 Done :)');
disp('Data has been channel/frame rejected and all parameter files have been generated.');
disp('Run ICA (AMICA or CUDAICA) by submitting .qsub or .sc file and then proceed to step 3 to fit dipoles from ICs.');
