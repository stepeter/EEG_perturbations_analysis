function [EEG,ICA_STRUCT]=interpolEEG(EEG_in,ICA_STRUCT_in,useExts,interpolate)
    %Interpolates missing channels in EEG while also reformatting
    %the ICA_STRUCT to reflect this interpolation
    
%    %Remove externals if specified 
%    if useExts==0
%        EEG_in=pop_select(EEG_in,'nochannel',(EEG_in.nbchan-7):EEG_in.nbchan);
%    end
   
   %Interpolate bad channels as listed in ICA_STRUCT_in
   if interpolate==1
       if isfield(ICA_STRUCT_in,'interpChans')
           targetChannels=ICA_STRUCT_in.interpChans;
           EEG = interpolateChannels(EEG_in, targetChannels); %default is to use all other channels for interp
           disp(['Interpolating ' num2str(length(targetChannels)) ' channels']);
           ICA_STRUCT=ICA_STRUCT_in;
       else
           targetChannels=setdiff(1:EEG_in.nbchan,ICA_STRUCT_in.good_chans);
           EEG = interpolateChannels(EEG_in, targetChannels); %default is to use all other channels for interp
           disp(['Interpolating ' num2str(length(targetChannels)) ' channels']); 
           
           %Re-format ICA_STRUCT (channel rejection history is still stored)
           ICA_STRUCT=ICA_STRUCT_in;
           ICA_STRUCT.good_chans=1:EEG.nbchan;
           ICA_STRUCT.good_cap=1:EEG.nbchan;
           ICA_STRUCT.interpChans=targetChannels; %save channels that were interpolated
       end
   else
       disp('No interpolation performed');
   end
