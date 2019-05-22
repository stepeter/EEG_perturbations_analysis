function EEG=reformatEventStruct(EEG,evSuffix)
%Adds event suffix to each event and changes type name to something more
%meaningful

for j=1:length(EEG.event)
    if regexp(EEG.event(1,j).type, regexptranslate('wildcard','* pressed'))
        %If it's a keypress event, change event type based on original type
        switch EEG.event(1,j).type(1) 
            case 'P'
                mainEv='M_on_CCW'; %counter-clockwise start
            case 'J'
                mainEv='M_off_CCW'; %counter-clockwise end
            case 'O'
                mainEv='M_on_CW'; %clockwise start
            case 'L'
                mainEv='M_off_CW'; %clockwise end
            otherwise
                disp('Other keypress detected! No action performed for this!');
                mainEv='other';
        end
        EEG.event(1,j).type=[mainEv '_' evSuffix];
    end
end
