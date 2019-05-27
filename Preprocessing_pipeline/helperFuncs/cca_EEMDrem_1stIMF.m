function [EEG,IMF1]=cca_EEMDrem_1stIMF(EEG)

EEG_orig=EEG;

%Run EEMD (can take awhile and can be a LOT of memory)
nmodes=1; %only need 1st mode 
[IMFs_final,EEG]=runEEMDlowmem(EEG,nmodes);

%Iterate across each IMF from EEMD
for i=1
    EEG_temp=EEG_orig;
    EEG_temp.data=squeeze(IMFs_final(:,i+1,:));
    
    %For first IMF (high frequencies): run CCA and
    %remove all components with higher IQR than median IQR of channel data
    %(just remove the big high frequency CCA components)
    [EEG_temp]=autoCCA_remEMG(EEG_temp); %runs pseudoinverse.m because low rank
    
    IMFs_final(1:size(EEG_temp.data,1),i+1,:)=EEG_temp.data;
end

IMF1=squeeze(IMFs_final(:,i+1,:));
%Sum up IMFs to get (ideally) cleaned data
EEG.data=squeeze(sum(IMFs_final(:,2:(nmodes+2),:),2)); %squeeze(sum(IMFs_final(:,2:(nmodes+1),:),2));
clear IMFs_final; %this can be a big variable, so remove it ASAP

% [EEG]=runAutoCCA_motion(EEG,badChanInds);


%Remove bad channels
% EEG=pop_select(EEG,'nochannel',badChanInds);
% pop_eegplot(EEG,1);
% figure; pop_spectopo(EEG, 1, [0 EEG.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
end



function [IMFs_final,EEG]=runEEMDlowmem(EEG,nmodes)

IMFs_final=zeros(EEG.nbchan,nmodes+2,size(EEG.data,2));

% tic
for i=1:EEG.nbchan
    disp(['Channel #' num2str(i)]);
    IMFs_final(i,:,:) = mod_pop_runemd(EEG, 'emdtype','EEMD','norm',0,'nmodes',nmodes,'ensemblenum',2,'noiseassist',0.1,'wsize',100,'ssize',20,'defaultnoisechannel',4,'dataset',1,'chanind',i);
end
EEG.imfnumber=nmodes;
EEG.eemdensemblenumber=2;
EEG.eemdassistednoise=0.1;
EEG.emdchansind=EEG.nbchan;
% toc
end

function [EEG]=autoCCA_remEMG(EEG)
%Run CCA and put it in ICA spot
delay=1;
[W,r] = bsscca_local(EEG.data,delay);

EEG.icaweights=[]; EEG.icasphere=[]; EEG.icawinv=[]; EEG.icaact=[];
EEG.icaweights=real(W);
EEG.icasphere=eye(EEG.nbchan);
EEG=eeg_checkset(EEG,'ica');

%Only remove components with higher IQR than median channel IQR
chanIQRthresh=median(iqr(EEG.data,2));
compIQRs=iqr(EEG.icaact,2);
inds=find(compIQRs>chanIQRthresh);
EEG=pop_subcomp(EEG,inds);
end

function [W,r] = bsscca_local(X,delay)
% bsscca() - Blind Source Separation through Canonical Correlation Analysis
% From AAR toolbox (EEGLAB extension)
%
% Usage:
%   >> [W,r] = bsscca(X,delay)
%
% Inputs:
%   X     - data matrix (dxN, data channels are rowwise)
%   delay - delay at which the autocorrelation of the sources will be
%           maximized (def: 1)
%
% Output:
%   W     - separation matrix
%   r     - autocorrelation of the estimated sources at the given delay
%
% See also:
%   BSSCCA_IFC, AUTOBSS

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 2, delay = 1; end
if nargin < 1, 
    help cca;
    return;
end


[d,T] = size(X);

% correlation matrices
Y = X(:,delay+1:end);
X = X(:,1:end-delay);
Cyy = (1/T)*Y*Y';
Cxx = (1/T)*X*X';
Cxy = (1/T)*X*Y';
Cyx = Cxy';


% calculate W
if rank(double(X'))<size(X,1)
    invCyy = pseudoinverse(Cyy);
    [W,r] = eig(pseudoinverse(Cxx)*Cxy*invCyy*Cyx);
else
    invCyy = pinv(Cyy);
    [W,r] = eig(pinv(Cxx)*Cxy*invCyy*Cyx);
end
r = sqrt(abs(real(r)));
[r,I] = sort(diag(r),'descend');
W = W(:,I)';
end

