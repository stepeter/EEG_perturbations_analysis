function [EEG EEG_VR Chan2upd ChanInfo] = rejPowerBase(EEG, ICA_STRUCT, Ctask)

EEG_VR = EEG;

%% re-referencing new dataset to average
EEG_VR = pop_reref( EEG_VR, []);

if isempty(ICA_STRUCT.rej_frame_idx) == 0
    
    %% defining sectors to exclude from the alternative EEG dataset
    sep_Frame = find(diff(ICA_STRUCT.rej_frame_idx) > 1)+1;
    
    for fr = 1:length(sep_Frame)
        if fr == 1
            sec(1,1) = ICA_STRUCT.rej_frame_idx(1); sec(1,2) = ICA_STRUCT.rej_frame_idx(sep_Frame(1)-1);
        elseif fr > 1
            sec(fr,1) = ICA_STRUCT.rej_frame_idx(sep_Frame(fr-1)); sec(fr,2) = ICA_STRUCT.rej_frame_idx(sep_Frame(fr)-1);
        end
    end
    
    % retrieving last period
    sec(fr+1,1) = ICA_STRUCT.rej_frame_idx(sep_Frame(fr)); sec(fr+1,2) = ICA_STRUCT.rej_frame_idx(end);
    
    FrRej = [];
    for vfg = 1:size(sec,1)
        FrRej = [FrRej num2str(sec(vfg,1)) ' ' num2str(sec(vfg,2)) ';'];
    end
    
    EEG_VR = pop_select(EEG_VR,'nopoint',str2num(FrRej));
end

figure; pop_spectopo(EEG_VR, 1, [0      EEG_VR.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');

conty = 1; Nrej = [];
while conty == 1
    display(['Current K value = ' num2str(EEG_VR.pnts/(EEG_VR.nbchan^2))])
    pause
    
    ImpName2 = input('Select channel range to exclude: ');
    EEG_VR = pop_select(EEG_VR,'nochannel',ImpName2);
    close all
    
    Vfig = figure;
    pop_spectopo(EEG_VR, 1, [0      EEG_VR.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
    pause
    
    conty = input('exclude more channels? (1 = yes, 0 - no): ');
end

%% saving figure of power spectrum
saveas(Vfig,['PowerSpectrum_' Ctask],'fig')

% %% saving cleaned dataset (with no bad frames)
% pop_saveset(EEG_VR,'filename', [Tsk(TaskType).name '_PwrAnalysis.set'],...
%     'filepath',EEGsets_outpath,'savemode','twofiles');

% getting channel labels (numbers) from real EEG dataset
for fg = 1:size(EEG.chanlocs,2)
    getF = EEG.chanlocs(1,fg);
    ChanInfo{fg,1} = {cellstr(char(getF.labels)), getF.urchan};
    Chan2comp(fg,1) = getF.urchan;clear getF
end

% getting channel labels (numbers) from cleaned EEG dataset
for fg = 1:size(EEG_VR.chanlocs,2)
    getF = EEG_VR.chanlocs(1,fg);
    Chan2comp_VR(fg,1) = getF.urchan;clear getF
end

%% finding channels to exclude from real dataset and excluding them
ChEx = find(ismember(Chan2comp,Chan2comp_VR) == 0); %% finds which channels from the adapted dataset are not present in the real dataset
EEG = pop_select(EEG,'nochannel',ChEx); clear Chan2*; clear fg;

% getting channel labels (numbers) from real EEG dataset
ChanInfo = [];
for fg = 1:size(EEG.chanlocs,2)
    getF = EEG.chanlocs(1,fg);
    ChanInfo{fg,1} = {cellstr(char(getF.labels)), getF.urchan};
    Chan2upd(1,fg) = getF.urchan;clear getF
end
