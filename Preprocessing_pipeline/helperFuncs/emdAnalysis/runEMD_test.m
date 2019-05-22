addpath(genpath('/media/stepeter/Local_Data/VR_rotations/Code'));

%Run EMD
tic
EEG = pop_runemd(EEG, 'emdtype','EEMD','norm',0,'nmodes',10,'ensemblenum',2,'noiseassist',0.1,'wsize',100,'ssize',20,'defaultnoisechannel',4,'dataset',1,'chanind',[1:128]);%[1 19 23 54 58 81 85 115 119] );
toc
EEG = eeg_checkset( EEG );

imfsRem=1; %[1 2]; %
%Remove IMF's
for i=1:length(imfsRem)
    EEG.data(EEG.emdchansind,:)=EEG.data(EEG.emdchansind,:)-squeeze(EEG.IMFs(:,imfsRem(i)+1,:));
end

%Plot EEG
pop_eegplot( EEG, 1, 1, 1);
%revert to original data: EEG.data=EEG_orig.data;

%Plot EMD results
%pop_eegplotemd( EEG, 2, 2, 1);
pop_spectopo(EEG);