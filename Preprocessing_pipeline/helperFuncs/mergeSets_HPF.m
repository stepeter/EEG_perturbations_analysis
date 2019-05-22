function [EEG]=mergeSets_HPF(EEG_processed,group,EEGsets_outpath)
%Set HPF and merge datasets

% cd(EEGsets_outpath);

for i=1:length(EEG_processed)
    %Load in the dataset
    EEG=EEG_processed(i);
    
    %Scan over all events and add suffix based on condition
    evSuffix=MISM_condAbbrevs(EEG.setname);
    EEG=reformatEventStruct(EEG,evSuffix);
    EEG_new(i)=EEG;
end

%Now merge the sets
EEG = [];
%EEG=mergeEEGsetsMISMCombs(EEG_new,group); %,set_files);
EEG = pop_mergeset(EEG_new, 1:length(EEG_new));
EEG.setname = ['Merge_' group];
EEG = eeg_checkset(EEG); % always checkset

%Add filter size to pathname
% pop_saveset(EEG, 'filename', ['Merge_' group], 'filepath', EEGsets_outpath,'savemode','twofiles');
% 
% %Plot power spectrum now that it is merged
% Vfig=figure;
% pop_spectopo(EEG, 1, [0 EEG.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
% saveas(Vfig,[EEGsets_outpath filesep 'PowerSpectrum' group],'fig');
% close all;
