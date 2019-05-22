function [cleanData bad_frames bad_frames_marked_zero index bad_frame_all_chans]= eeg_badframes(EEG,iqr_average)
% [cleanData bad_frames bad_frames_marked_zero index bad_frame_all_chans]=
% eeg_badframes(EEG,iqr_average)
%
%input EEG can be an EEG set OR a data matrix

%modifications
%JTG added index and bad_frames_all_chans to output 5/18/09 
%   index is like a normalized power
%   bad_frames_all_chans = mean(index')
%JTG removed print to screen
%JTG 10/02/09 first input can be EEG set OR data matrix

if nargin<2
    iqr_average=0.8;
end; 

if isstruct(EEG)
    fprintf('calculating power, median power and power iqr...\n');
    power = (EEG.data(:,:) .^2)'; % for each sample point and each channel seperately power is calculated
    nbchan = EEG.nbchan;
    cleanData = EEG.data(:,:);
else
    fprintf('calculating power, median power and power iqr...\n');
    power = (EEG.^2)'; % for each sample point and each channel seperately power is calculated
    nbchan = size(EEG,1);
    cleanData = EEG;
end
pow_median=median(power); % calcuate the median of power for each channel seperately
pow_iqr = iqr(power);

fprintf('calculating power iqr distance to its median for each frame...\n');
h = waitbar(0,'please wait...');
for i=1:nbchan
    index(:,i) = ((power(:,i) - pow_median(i))/pow_iqr(i)); 
    % kind of a standardization (z-transform) but taking power-median and
    % power inter-quartile distance instead of mean and sd
    waitbar(i/nbchan);
end;
close(h);
% now the mean of the index value (standardized power) for all channels is
% computed
bad_frame_all_chans = mean(index');
%fprintf('%d percent of the frames are bad.\n',100*length(find(bad_frame_all_chans>iqr_average))/length(bad_frame_all_chans));
bad_frames = find(bad_frame_all_chans>iqr_average); % if mean is over a certain value it is defined as bad
bad_frames_marked_zero = ones(1,size(bad_frame_all_chans,2));
bad_frames_marked_zero(bad_frames) = 0;


cleanData(:,bad_frames) = [];

return

% rey's way
pow=sqrt(sum(data.^2));
pow_median=median(pow);
pow_iqr=iqr(pow);
bad_frames = find(pow > pow_median+5*pow_iqr); % was find(pow > pow_median+5*pow_iqr) before
bad_frames_marked_zero = ones(1,size(data,2));
bad_frames_marked_zero(bad_frames) = 0;