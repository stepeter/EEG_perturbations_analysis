function [return_array idx] = get_EEG_event_array(EEG,type,field);

%function [return_array idx] = get_EEG_event_array(EEG,type,field);
%
%returns an array of the values in EEG.event."field" for any event that 
%matches input "type"
%type can be a string or an integer

%J Gwin 05/21/09
return_array = [];
idx = [];

%if type is a cell then convert
if iscell(type)
    type = type{1};
end

if isnumeric(type)
    n = 0;
    for i = 1:length(EEG.event)
        if EEG.event(i).type == type
            n = n + 1;
            eval(['return_array{n} = EEG.event(i).' field ';']);
            idx(n) = i;
        end
    end
elseif ischar(type)
    n = 0;
    for i = 1:length(EEG.event)
        if strcmp(EEG.event(i).type,type)
            n = n + 1;
            eval(['return_array{n} = EEG.event(i).' field ';']);
            idx(n) = i;
        end
    end
else
    error('Type should be a string or numeric');
end
