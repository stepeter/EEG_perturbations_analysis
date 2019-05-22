function [filenames fullfilenames] = dir_list(filepath,filetype)
%function filenames = dir_list(filepath,filetype)
%returns a cell array of all filenames and sub-directory names
%in a given folder
%
%inputs:
%   filepath - directory as a string
%   filetype - optional input variable if defined function only returns
%              files of thies type otherwise returns all files and 
%              sub-folders. Filetype should be a string (ex '.xls' or 'xls')
%
%
%Created by JTG 10/12/06

%Modifications

%Copyright Simbex 2006
filenames = [];

if nargin == 1
    struct = dir(filepath);
    n = 0;
    for i = 1:length(struct)
        if ~strcmp(struct(i).name,'.') & ~strcmp(struct(i).name, '..')
            n = n + 1;
            filenames{n} = struct(i).name;
            fullfilenames{n} = fullfile(filepath,struct(i).name);
        end
    end
else
    struct = dir(filepath);
    n = 0;
    for i = 1:length(struct)
        if struct(i).isdir == 0
            if strcmp(struct(i).name(end-(length(filetype)-1):end), filetype)
                n = n + 1;
                filenames{n} = struct(i).name;
                fullfilenames{n} = fullfile(filepath,struct(i).name);
            end
        end
    end 
end