if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_connectivity/Code/EEG/'));
    startUpEEGLAB('close');
end
close all;


for subjnum=[1:13 15:17 19:20 22:33]
    %Set paths
    RootData='/media/stepeter/Local_Data/VR_connectivity/Data/';  
    subjcode = ['WMISM_' num2str(subjnum)];
    group='all'; %'mismatches'; %'pull';
    subj_basepath = [RootData subjcode filesep]; % Top level subject folder
    EEGsets_outpath = [subj_basepath 'EEG_sets']; % Where you want to save your .set files
%     ICA_path = [subj_basepath 'ICA_Stuff'];
    ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
    
    %Load set and ICA_STRUCT
    cd(subj_basepath)
    global ALLEEG EEG CURRENTSET ALLCOM;
    EEG=pop_loadset('filename', ['Merge_' group '.set'], 'filepath', ['/media/stepeter/Local_Data/VR_rotationsDenoise/Data/' subjcode filesep 'EEG_sets']);    
    load(['/media/stepeter/Local_Data/VR_rotationsDenoise/Data/' subjcode filesep 'ICA_Stuff' filesep 'files_AMICA' group filesep group '_ICA_PREP_1.mat']);
    EEG_orig=EEG; EEG_orig=pop_select(EEG_orig,'channel',1:128);
    
    %Select training & pre/post trials only
    segmentBounds='boundary';
    EEGin=EEG;
    
    %Pull out training and pre/post trials
    if ischar(segmentBounds) && strcmp(segmentBounds,'boundary')
        %Get boundary events
        boundary_latencies=[];
        for i=1:length(EEGin.event)
            if strcmp(EEGin.event(i).type,'boundary')
                boundary_latencies=[boundary_latencies EEGin.event(i).latency];
            end
        end
        %     %Special condition for subject data
        %     if strcmp(EEGin.subject,'WMISM_1')
        %         boundary_latencies(end-1)=[];
        %     end
        boundary_latencies2=[1 boundary_latencies EEGin.pnts];
    else
        boundary_latencies2=segmentBounds;
    end
    EEG=pop_select(EEGin,'point',boundary_latencies2([1 5]));
    
    %Remove bad channels    
    EEG=pop_select(EEG,'channel',ICA_STRUCT.good_chans); %remove bad channels
    
    %Denoising Methods:
    %1) ASR to remove ungrounding events (use st. dev. of 20)
    EEGtemp=clean_asr(EEG,20); EEG.data=EEGtemp.data; EEGtemp=[];
    
    %2) Run EEMD and CCA to remove large components in 1st IMF
    [EEG,IMF1]=cca_EEMDrem_1stIMF(EEG);
    EEGtemp=EEG; EEGtemp.data=IMF1; %plot the first IMF
    Vfig=figure; pop_spectopo(EEGtemp, 1, [0 EEGtemp.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 100],'electrodes','off');
    saveas(Vfig,[EEGsets_outpath filesep 'IMF1_spec'],'fig');
    saveas(Vfig,[EEGsets_outpath filesep 'IMF1_spec'],'png');
    
%     %3) Use EOG externals to remove eye activity using ICA
%     [wts,sph]=binica(EEG.data,'extended',1); %run ICA
% %     %Load in old results instead
% %     wtsFile=dir([subj_basepath '*.wts']);
% %     wts=floatread([subj_basepath wtsFile.name], [1 Inf]);
% %     wts=reshape(wts,sqrt(length(wts)),sqrt(length(wts)));
% %     sphFile=dir([subj_basepath '*.sph']);
% %     sph=floatread([subj_basepath sphFile.name], [1 Inf]);
% %     sph=reshape(sph,sqrt(length(sph)),sqrt(length(sph)));
%     
%     EEG.icaweights=wts; EEG.icasphere=sph; EEG=eeg_checkset(EEG,'ica');
%     [EEG,eye_artifacts]=interface_ADJ(EEG,'test.txt'); %run Adjust plugin to find eye artifacts
%     %Also remove each IC that has maximum weight for neck externals
%     [~,maxEMGIC1]=max(abs(EEG.icawinv(EEG.nbchan-4,:)));
%     [~,maxEMGIC2]=max(abs(EEG.icawinv(EEG.nbchan-3,:)));
%     %Also remove each IC that has maximum weight for eye externals (likely
%     %removed by Adjust already)
%     [~,maxEOGIC1]=max(abs(EEG.icawinv(EEG.nbchan-2,:)));
%     [~,maxEOGIC2]=max(abs(EEG.icawinv(EEG.nbchan-1,:)));
%     [~,maxEOGIC3]=max(abs(EEG.icawinv(EEG.nbchan,:)));
%     
%     final_artifacts=unique([eye_artifacts maxEMGIC1 maxEMGIC2 maxEOGIC1 maxEOGIC2 maxEOGIC3]);
%     EEG=pop_subcomp(EEG,final_artifacts); %remove eye components
    EEG.icaweights=[]; EEG.icasphere=[]; EEG.icawinv=[]; EEG.icaact=[]; EEG.icachansind=[]; %clear out ICA info
    
    %Remove externals (last 5 channels)
    EEG=pop_select(EEG,'nochannel',(EEG.nbchan-4):EEG.nbchan);
    
    %Re-reference and interpolate before saving
    EEG=pop_reref(EEG,[]);
    EEG=pop_interp(EEG,EEG_orig.chanlocs,'spherical');
    EEG.setname=ICA_STRUCT.associated_set;
    
    pop_saveset(EEG,'filename', ['Merge_' group '_CAR_v2.set'], 'filepath', EEGsets_outpath);
    %Plot the output
%     pop_eegplot(EEG,1);
    Vfig=figure; pop_spectopo(EEG, 1, [0 EEG.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 100],'electrodes','off');
    saveas(Vfig,[EEGsets_outpath filesep 'Powerspec_CAR_v2'],'fig');
    saveas(Vfig,[EEGsets_outpath filesep 'Powerspec_CAR_v2'],'png');
    close all;
end