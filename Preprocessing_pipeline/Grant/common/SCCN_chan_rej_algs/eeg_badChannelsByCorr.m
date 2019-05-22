function badChans = eeg_badChannelsByCorr(EEG, percentTimeBad)
% find bad channels using correlation with nearby channels and mark them
% for rejection if they have gone bad more than 'percentTimeBad' of
% experiment time
% Usage:
%   badChans = eeg_badChannelsByCorr(EEG, percentTimeBad)

if nargin < 2
    percentTimeBad = 0.01;
end;

[badchanindts,badchanlogts,c_cellts,IND,maxcts]=detectnoisychannelts(EEG);

percentSecondsBad = sum(badchanlogts') / size(badchanlogts,2);

badChans = find(percentSecondsBad>percentTimeBad);
