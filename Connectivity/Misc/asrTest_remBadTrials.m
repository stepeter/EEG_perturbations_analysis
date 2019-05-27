%For non-ASR data, remove bad trials with abnormally high amplitude

loadpath='/usr/local/VR_connectivity/Data/regularGroupConn_ASRtest/conns/orig/';
outpath='/usr/local/VR_connectivity/Data/regularGroupConn_ASRtest/conns/';
epochTimes=[-1 2];
EpPre = epochTimes(1,1)+(-0.6);
EpPost = epochTimes(1,2)+0.6;
%%
A=dir([loadpath '*.set']);

for i=1:length(A)
    EEG=pop_loadset('filepath',loadpath,'filename',A(i).name);
    %Get rid of really noisy epochs (due to ungrounding)
    EEGtemp=EEG; EEGtemp.data((end-7):end,:)=0; EEGtemp.icaact((end-7):end,:)=0;
    if i==1
        [~,inds]=pop_eegthresh(EEGtemp,1,1:EEG.nbchan,-500,500,EpPre,EpPost,0,1);
    else
        [~,inds]=pop_eegthresh(EEGtemp,1,1:EEG.nbchan,-2000,2000,EpPre,EpPost,0,1);
    end
    EEG2=pop_select(EEG,'notrial',inds);
    pop_saveset(EEG2,'filepath',outpath,'filename',A(i).name);
end