function EEG=removeExtraKeyboardEvents(EEG)
%%Remove Extra EEG events from releasing keyboard keys

%Find event indices to remove
inds=[];
for i=1:length(EEG.event)
    if regexp(EEG.event(1,i).type, regexptranslate('wildcard','* released')) %removes all events with '* released' label
        inds=[inds, i];
    end
end

%Now remove those events
EEG=pop_editeventvals(EEG,'delete',inds);
EEG = eeg_checkset(EEG,'eventconsistency');
end