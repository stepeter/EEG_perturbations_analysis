% clear; close all; clc;

%Set parameters
subjsToAnalyze=[1:13 15:17 19:20 22:33]; %for batch processing
group='all'; %grouping to look at
numChans=136;

%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_rotationsDenoise/Code/'));
    startUpEEGLAB('close');
end

%Define folder where data is stored (.xdf files)
pathname='/media/stepeter/Local_Data/VR_rotationsDenoise/Data/';
cd(pathname);

for i=subjsToAnalyze
    %Load in the data
    remExts=0; %remove externals (0 - keep exts. 3-7, which were not mastoid or EKG)
    [EEG_raw,subj_basepath]=EEG_load_beforePREP(i,pathname,numChans,'common',group,remExts);
    
    %Set up path to save .set files
    EEGsets_outpath = [subj_basepath 'EEG_sets'];
    if ~exist(EEGsets_outpath)
        mkdir(EEGsets_outpath)
    end
    
    %Downsample, high-pass filter, and run cleanline
    for d=1:length(EEG_raw)
        %Downsample to 256Hz
        setname=EEG_raw(d).setname;
        EEG_raw(d) = pop_resample(EEG_raw(d), 256);
        EEG_raw(d).setname=setname; %revert to regular setname
        
        %High-pass filter at 1 Hz
        filter_hi=0; filter_lo=1;
        EEG_raw(d) = pop_eegfiltnew(EEG_raw(d),filter_lo,filter_hi);
        
        %Common median reference (not on externals)
        refMedian=median(EEG_raw(d).data(1:128,:));
        EEG_raw(d).data=EEG_raw(d).data-ones(EEG_raw(d).nbchan,1)*refMedian;
        
        %Run cleanline
        EEG_processed(d) = pop_cleanline(EEG_raw(d), 'bandwidth', 2,'chanlist', [1:EEG_raw(d).nbchan], 'computepower', 0, 'linefreqs', [60 120 180 240],...
            'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100,...
            'verb', 1, 'winsize', 4, 'winstep', 4);   
    end
    
    %Merge the sets
    [EEG]=mergeSets_HPF(EEG_processed,group,EEGsets_outpath); 
end

disp('Step 1 Done :)');
disp('Data has been formatted, high-pass filtered, and merged.');
disp('Proceed to step 2 to remove bad channels/frames.');
