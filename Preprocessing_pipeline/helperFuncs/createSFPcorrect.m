function [ELocs_file_eeglab]=createSFPcorrect(ELocs_file,numChans,ElecFolder)
%Creates the sfp file with channel locations needed for eeglab
%numChans: number of channels (not including externals)

% Use Head models 
if ~exist([ELocs_file(1:end-4) '_eeglabformat_final.sfp'])
    ELocs_rawfile=ELocs_file;
    ELocs_file=[ELocs_rawfile(1:(length(ELocs_rawfile)-4)) '_eeglabformat.sfp'];
    standardcapelocs=['/home/stepeter/HNL_Cluster/share/data3/stepeter/Steve_Code/Biosemi/cap128.sph'];
    %standardcapelocs='/home/stepeter/HNL_Cluster/share/data3/stepeter/Grant/common/Biosemi_original_files/cap256.sph'; %use template used in zebris digitization
    collectionsys='Biosemi';
    [ELocs_file_eeglab]=hnl_convertZebris_WMISM(ELocs_rawfile,ELocs_file, standardcapelocs,collectionsys,numChans,ElecFolder);
else
    ELocs_file_eeglab = [ELocs_file(1:end-4) '_eeglabformat_final.sfp'];
end
