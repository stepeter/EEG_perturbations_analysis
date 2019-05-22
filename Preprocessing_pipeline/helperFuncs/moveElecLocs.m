function EEG=moveElecLocs(EEG,fromChans,toChans)
%%Get chanloc information from fromChans and add to toChans without
%%changing toChans labels
if length(fromChans)~=length(toChans)
    error('fromChans and toChans need to be the same size!');
end
%Get fromChans channel locations
temp=EEG.chanlocs(1,fromChans);

%Get toChan labels and replace toChans info with updated values
q=0;
for p=toChans
    q=q+1;
    temp(1,q).labels=EEG.chanlocs(1,p).labels;
    temp(1,q).urchan=EEG.chanlocs(1,p).urchan;
    temp(1,q).type=EEG.chanlocs(1,p).type;
    EEG.chanlocs(1,p)=temp(1,q);
end
end