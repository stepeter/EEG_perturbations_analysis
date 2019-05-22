function sets2Load=mergeEEGsetsMISMCombs(group,subjnum) %,set_files)
%%Takes in all EEG sets from MISM experiment and only merges a subset of
%%those based on the group parameter ('mismatches' or 'motor_learning')

switch group
    case 'motor_learning'
        sets2Load={'Pre','Train1','Train2','Train3','Post'};
    case 'mismatches'
        sets2Load={'Stand_Viz','Walk_Viz'}; %,'Balance_noMism','Balance_Mism'};
    case 'pull'
        sets2Load={'Stand_Pull','Walk_Pull'};
    case 'all'
        if subjnum==1
            sets2Load={'Stand_Pull','Walk_Pull','Stand_Viz','Walk_Viz','Pre','Train1','Train2','Train3','Train3_2','Post'};
        else
            sets2Load={'Stand_Pull','Walk_Pull','Stand_Viz','Walk_Viz','Pre','Train1','Train2','Train3','Post'};
        end
    otherwise
        error('Group name not recognized!');
end
