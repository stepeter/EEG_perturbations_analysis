function [finalelocs]=hnl_convertZebris(ELocs_rawfile,ELocs_file,standardcapelocs,collectionsys)
%Convert Zebris digitizer data into EEGlab format 

%Inputs:    
%   ELocs_rawfile: .sfp file from Zebris digitzer output
%   ELocs_file: final file to be imported into dataset
%   standardcapelocs: standard file from Biosemi cap
%   collectionsys: EEG collection system used, 'DataRiver' or 'Biosemi'

%A Sipp 05-24-11
 
%Determine correct electrode labeling/naming:
% [whatdigitized,ok] = listdlg('PromptString','What electrode locations were digitized?',...
%             'SelectionMode','single','ListString',{'standard 256 + 8 externals'; 'standard 256'; ...
%             'Sample of 256'; 'Sample of 256 + 8 externals'}, 'InitialValue', 4,'ListSize',[400 100]);

        
whatdigitized=1;
ok=1;
        
        
        
%Get standard labels
stanlocs=fopen(standardcapelocs);
C = textscan(stanlocs, '%f %f %f %s');
fclose(stanlocs);
if strcmp(collectionsys,'Cognionics') == 0
    stand_labels=C{4}(1:256);
else
    stand_labels=C{4}(4:67);
end
if strcmp(collectionsys,'DataRiver') == 1
    stand_labels=[stand_labels' 'I1' 'I2' 'I3' 'I4' 'I5' 'I6' 'I7' 'I8'];
else
    stand_labels=[stand_labels' 'Ext1' 'Ext2' 'Ext3' 'Ext4' 'Ext5' 'Ext6' 'Ext7' 'Ext8'];
end

%Get sample labels
sampleelocs='A1,A17,A23,B1,B25,B30,C12,C21,C25,D4,D12,D27,E16,E20,E28,F2,F23,F28,G3,G18,G22,H14,H18,H22';
if whatdigitized == 3 || whatdigitized == 4
    [elocsample,ok] = listdlg('PromptString','What were the sample electrode loctions (only include 1-256)?',...
            'SelectionMode','single','ListString',{sampleelocs; 'other'}, 'InitialValue', 1,'ListSize',[800 100]);
    if elocsample == 1
        sampleelocs={'A1','A17','A23','B1','B25','B30','C12','C21','C25','D4','D12','D27','E16','E20','E28','F2','F23','F28','G3','G18','G22','H14','H18','H22'};
    else
        sampleelocs=input('Sample electrode locations: ', 's');
    end
end

%Create digitizer data label list that will correspond to EEG recording software
if whatdigitized == 1 %'standard 256 + 8 externals'
    dig_labels = stand_labels; %these are the channel names that will be in raw zebris output file if a full set is digitized
elseif whatdigitized == 2 %'standard 256 '
    dig_labels = stand_labels(1:256); %these are the channel names that will be in raw zebris output file if a full set is digitized
elseif whatdigitized == 3 % 'Sample of 256'
    dig_labels=sampleelocs;
elseif whatdigitized == 4 % 'Sample of 256 + 8 externals'
    dig_labels=[sampleelocs stand_labels(257:264)];
end

%%%% %LPA, NZ, RPA = 'fidt9','fidnz','fidt10' %%%
fiducial_labels={'LPA', 'NZ', 'RPA'};

%load digitizer data file
rawlocs=fopen(ELocs_rawfile);
rawC = textscan(rawlocs, '%s %f %f %f',3);
rawE = textscan(rawlocs, '%f %s %f %f %f');
fclose(rawlocs);

%Fiducal coordinates
fidchanlabel=rawC{1};
fidXcoord=rawC{2};
fidYcoord=rawC{3};
fidZcoord=rawC{4};

%Electrode coordinates
elocchanlabel=[rawE{2}];
elocXcoord=[rawE{3}];
elocYcoord=[rawE{4}];
elocZcoord=[rawE{5}];

%Convert to m from raw Zebris data (mm)
% elocXcoord=elocXcoord/1000;
% elocYcoord=elocYcoord/1000;
% elocZcoord=elocZcoord/1000;

%create seperate file for Fiducal and Electrode locations
%ELocs_file='/share/data/amysipp/Amy/CL2/Cognitive_Loading_CL2_eeglabformat2.sfp';


if ~exist([ELocs_file(1:(length(ELocs_file)-4)) '_final.sfp'])
    if ~isempty(ELocs_file)==1
        %Digitized electrode locations file
        fid = fopen(ELocs_file,'w');
        for n = 1:length(elocXcoord)
            fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\t%i\n',dig_labels{n},elocXcoord(n),elocYcoord(n),elocZcoord(n), n);
        end
        fclose(fid);
        disp(['Saved ' ELocs_file]);
        
        %Digitized fiducal locations file
        fid = fopen([ELocs_file(1:(length(ELocs_file)-4)) '.fid'],'w');
        for n = 1:length(fidXcoord)
            fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\t%i\n',fiducial_labels{n},fidXcoord(n),fidYcoord(n),fidZcoord(n),n);
        end
        fclose(fid);
        disp(['Saved ' ELocs_file(1:(length(ELocs_file)-4)) '.fid']);
    end
    
    plo=1;
    if plo ==1
        digHCS=[elocXcoord elocYcoord elocZcoord];
        % digHCS=[elocXcoord(1:31) elocYcoord(1:31) elocZcoord(1:31)];
        figure; set(gcf,'color','w');
        plot3(digHCS(:,1),digHCS(:,2),digHCS(:,3),'rx');
        hold on; plot3([0;.05],[0;0],[0;0],'r-','linewidth',2);
        hold on; plot3([0;0],[0;.05],[0;0],'g-','linewidth',2);
        hold on; plot3([0;0],[0;0],[0;.05],'b-','linewidth',2);
        title('Digitized Electrode Locations');
    end
    
    %Morph sample eloc set to full eloc set
    %don't need this if digitizing all electrodes:
    % if nose is just in x of y wrong, can use: with +X, -X, +Y, or -Y
    % as options
    %EEG=pop_chanedit(EEG, 'nosedir','+Y');
    % if ~isempty(dig_labels) && length(dig_labels) < 256
    %     copyfile([ELocs_file(1:(length(ELocs_file)-4)) '.fid'],[ELocs_file(1:(length(ELocs_file)-4)) '_final.fid'])
    %     copyfile(ELocs_file,[ELocs_file(1:(length(ELocs_file)-4)) '_final.sfp'])
    %     digHCS=[elocXcoord elocYcoord elocZcoord];
    %     MorphDigitizedSampleELocstoFullCap(standardcapelocs,[ELocs_file(1:(length(ELocs_file)-4)) '_final.sfp'],digHCS,dig_labels, 1,collectionsys)
    % end
    copyfile([ELocs_file(1:(length(ELocs_file)-4)) '.fid'],[ELocs_file(1:(length(ELocs_file)-4)) '_final.fid'])
    copyfile(ELocs_file,[ELocs_file(1:(length(ELocs_file)-4)) '_final.sfp'])
end
    
finalelocs=[ELocs_file(1:(length(ELocs_file)-4)) '_final.sfp'];
fprintf([finalelocs ' done \n']);
    
   



    