function [EEG]=trigChan2syncEvs(EEG)
%%Uses the trigger channel from xdf EEG set and converts it to events
%%before removing the channel. Then, you don't have to count channels
%%weirdly since the trigger channel is the first channel. Also removes the
%%16 AUX channels at the end.

%Find sync times (for a
diffTrig=diff(EEG.data(1,:));
syncRiseEvs=find(diffTrig>18 & diffTrig<50);
syncFallEvs=find(diffTrig<-18 & diffTrig>-50);

%Define terms
oldEvLen=length(EEG.event);
newEvLen=oldEvLen+length(syncRiseEvs)+length(syncFallEvs);

%Concatenate matrices of rising and falling edges
syncs=[syncRiseEvs syncFallEvs];
types=[zeros(1,length(syncRiseEvs)) ones(1,length(syncFallEvs))]+1;
typeNames={'Sync Rising','Sync Falling'};

% Store old trigger events into EEG.trigger
EEG.trigger=EEG.event;

%Add in each new event
q=0;
for i=oldEvLen+1:newEvLen
    q=q+1; %counter variable
    EEG.event(i).type=typeNames{types(q)};
    EEG.event(i).latency=(EEG.times(syncs(q))/1000*EEG.srate)+1; %converts from time (ms) to pnts
    EEG.event(i).duration=1;
end
EEG = eeg_checkset(EEG,'eventconsistency');

% %Remove data before first sync and after last sync
% syncsSorted=sort(syncs,2,'ascend');
% firstCut=EEG.times(syncsSorted(1))/1000-0.25; %quarter second before sync
% lastCut=EEG.times(syncsSorted(end))/1000+0.25; %quarter second after sync
% if (lastCut<(EEG.times(end)/1000))
%     EEG = pop_select(EEG,'notime',[lastCut EEG.times(end)/1000]); %Cut off end
% end
% EEG = pop_select(EEG,'notime',[EEG.times(1)/1000 firstCut]); %Cut off beginning

%Remove trigger channel (save as EEG.triggerchan)
EEG.triggerchan=EEG.data(1,:);
EEG=pop_select(EEG,'nochannel',1);
EEG=eeg_checkset(EEG);

%Remove last 16 AUX channels too
EEG=pop_select(EEG,'nochannel',(EEG.nbchan-15):EEG.nbchan);
EEG=eeg_checkset(EEG);