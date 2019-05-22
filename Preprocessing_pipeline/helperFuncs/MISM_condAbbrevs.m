function [abbrev]=MISM_condAbbrevs(cond)
%%Takes condition and changes it to specified abbreviation (case specific!)
switch cond
    case 'Stand_Pull'
        abbrev='SPL';
    case 'Walk_Pull'
        abbrev='WPL';
    case 'Stand_Viz'
        abbrev='SVZ';
    case 'Walk_Viz'
        abbrev='WVZ';
    case 'Pre'
        abbrev='PRE';
    case 'Post'
        abbrev='TPO';
    case 'Train1'
        abbrev='1T';
    case 'Train2'
        abbrev='2T';
    case 'Train3'
        abbrev='3T';
    case 'Train3_2'
        abbrev='3T';
    otherwise
        error('Unexpected filename!');
end
