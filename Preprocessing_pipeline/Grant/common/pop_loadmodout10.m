% pop_loadmodout() - import AMICA results
%
% Usage:
%   >> OUTEEG = pop_loadmodout( EEG, filepath );
%
% Inputs:
%   EEG            - input dataset
%   filepath       - file path for loadmodout function
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme & Jason palmer, SCCN, INC, UCSD

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: pop_loadmodout.m,v $

function [EEG, command] = pop_loadmodout10(EEG, filepath); 

    if nargin < 1
        help pop_loadmodout;
        return;
    end;
    
    command = '';
    if nargin < 2
        % ask user
        filepath = uigetdir('*.*', 'Choose a AMICA output folder -- pop_loadmodout'); 
        if filepath == 0 return; end;
    end;
    
    % import data
    % -----------
    mod = loadmodout10([ filepath filesep ]);

    EEG.icaweights = mod.W;
    EEG.icasphere  = mod.S;
    EEG.icawinv    = mod.A;
    EEG.icaact     = [];
    disp('AMICA weights imported successfully');
    
