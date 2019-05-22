%mod_std_precomp() JGwin modification to std_precomp
%   added feature so that if erspparms value corresponding to 'timewarp'
%   field is 'subject tw matrix' then I use a tw matrix that is stored in
%   ALLEEG.timewarp.latencies and similarly for 'timewarpms' if value is
%   'subject tw matrix' I use ALLEEG.timewarp.new_latencies
%
%std_precomp() - Precompute measures (ERP, spectrum, ERSP, ITC) for channels in a study. 
%                 If channels are interpolated before computing the measures, the updated 
%                 EEG datasets are also saved to disk. Called by pop_precomp(). Follow with 
%                 pop_plotstudy(). See Example below.
% Usage:    
% >> [ALLEEG,STUDY] = std_precomp(STUDY, ALLEEG, chanorcomp, 'key', 'val', ...);
%
% Required inputs:
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   chanorcomp   - ['components'|'channels'| or channel cell array] The string 
%                  'components' forces the program to precompute all selected 
%                  measures for components. The string 'channels' forces the 
%                  program to compute all measures for all channels.
%                  A channel cell array containing channel labels will precompute
%                  the selected measures. Note that the name of the channel is
%                  not case-sensitive.
% Optional inputs:
%  'erp'      - ['on'|'off'] pre-compute ERPs for each dataset.
%  'spec'     - ['on'|'off'] pre-compute spectrum for each dataset.
%               Use 'specparams' to set spectrum parameters.
%  'ersp'     - ['on'|'off'] pre-compute ERSP for each dataset.
%               Use 'erspparams' to set time/frequency parameters.
%  'itc'      - ['on'|'off'] pre-compute ITC for each dataset.
%               Use 'erspparams' to set time/frequency parameters.
%  'scalp'    - ['on'|'off'] pre-compute scalp maps for components.
%  'allcomps' - ['on'|'off'] compute ERSP/ITC for all components ('off'
%               only use pre-selected components in the pop_study interface).
%  'specparams' - [cell array] Parameters for the spectopo function are given as 
%              optional arguments:
%                   'freqrange' = [min max] frequency range to calculate. Changes x-axis limits {default: 
%                                1 Hz for the min and Nyquist (srate/2) for the max. If specified 
%                                power distribution maps are plotted, the highest mapped frequency 
%                                determines the max freq}. Note that it is better here to compute
%                                spectrum over a wide range of frequencies (it will then
%                                be possible to select another subrange for plotting).
%                   'freqfac'  = [integer] ntimes to oversample -> frequency resolution {default: 2}
%                   'nfft'     = [integer] length to zero-pad data to. Overwrites 'freqfac' above.
%                   'winsize'  = [integer] window size in data points {default: from data}
%                   'overlap'  = [integer] window overlap in data points {default: 0}
%                   'percent'  = [float 0 to 100] percent of the data to sample for computing the 
%                                spectra. Values < 100 speed up the computation. {default: 100}.
%                   'mapnorm'  = [float vector] If 'data' contain the activity of an independant 
%                                component, this parameter should contain its scalp map. In this case
%                                the spectrum amplitude will be scaled to component RMS scalp power.
%                                Useful for comparing component strengths {default: none}
%                   'rmdc'     = ['on'|'off'] 'on' -> remove DC {default: 'off'}  
%
%              Note that it is advised to compute spectrum 
%              over all frequencies since plotting function can always reduce
%              the range of plotted frequencies.
%  'erspparams' - [cell array] Optional arguments are 'cycles', 'freqrange',
%              'padratio', 'winsize', 'alpha' (see newtimef()). Note that it 
%              is adivised to select the largest frequency range and time window
%              as plotting function are capable of plotting subranges of
%              these. An important optional parameter that is
%                    'savetrials' = ['on'|'off'] save single-trials ERSP.
%                                   Requires a lot of disk space (dataset
%                                   space on disk times 10) but allow for
%                                   refined single-trial statistics.
%  'recompute' - ['on'|'off'] force recomputing ERP file even if it is 
%                already on disk.
%
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified by adding preprocessing 
%                  data as pointers to Matlab files that hold the pre-clustering component measures.
%   STUDY        - the input STUDY set with pre-clustering data added, for use by pop_clust() 
%
% Example:
%   >> [ALLEEG STUDY] = std_precomp(STUDY, ALLEEG,  { 'cz' 'oz' }, 'interpolate', 'on', 'erp', 'on', ...
%          'spec', 'on', 'ersp', 'on', 'erspparams', { 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });
%                          
%           % This prepares, channels 'cz' and 'oz' in the STUDY datasets.
%           % If a data channel is missing in one dataset, it will be
%           % interpolated (see eeg_interp()). The ERP, spectrum, ERSP, and 
%           % ITC for each dataset is then computed. 
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006, arno@sccn.ucsd.edu
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

% $Log: std_precomp.m,v $
% Revision 1.29  2008/11/13 02:45:52  arno
% better interpolation
%
% Revision 1.28  2008/11/13 02:45:23  arno
% revert version 1.26
%
% Revision 1.26  2008/11/13 02:34:54  arno
% nothing
%
% Revision 1.25  2008/11/12 23:09:25  arno
% Matlab 6.5 compatibility
%
% Revision 1.24  2008/04/19 21:02:25  arno
% fix problem concerning computation of channel indices
%
% Revision 1.23  2008/04/16 18:40:32  arno
% do not initialize changrp for component
% computation
%
% Revision 1.21  2008/02/15 16:51:36  arno
% simplify code for merging channel location files
%
% Revision 1.20  2007/12/09 00:40:15  arno
% recompute for topo
%
% Revision 1.19  2007/11/22 23:34:55  arno
% header
%
% Revision 1.18  2007/11/22 23:34:03  arno
% header
%
% Revision 1.17  2007/11/21 16:40:53  arno
% help msg
%
% Revision 1.16  2007/10/25 00:59:39  nima
% spectopo parameters described in help message.
%
% Revision 1.15  2007/09/11 10:51:16  arno
% precompute measures for components
%
% Revision 1.14  2007/04/06 22:09:44  arno
% recompute tag
%
% Revision 1.13  2007/04/05 23:17:55  arno
% guimode
%
% Revision 1.12  2007/04/05 23:13:39  arno
% *** empty log message ***
%
% Revision 1.11  2007/02/28 12:05:14  arno
% option to force recomputation
%
% Revision 1.6  2007/01/29 10:50:27  arno
% fix ERSP options
%
% Revision 1.4  2006/11/14 04:12:53  arno
% [Asame
%
% Revision 1.3  2006/11/14 03:59:25  arno
% debug ERSP check
%
% Revision 1.2  2006/11/14 03:53:18  arno
% Now checking file on disk
%
% Revision 1.1  2006/09/12 18:43:54  arno
% Initial revision
%

function [ STUDY, ALLEEG ] = mod_std_precomp_v10_2_5_5a(STUDY, ALLEEG, chanlist, varargin)
    
    if nargin < 2
        help std_precomp;
        return;
    end;
    
    if nargin == 2
        chanlist = 'channels'; % default to clustering the whole STUDY 
    end   
    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end

    g = finputcheck(varargin, { 'erp'         'string'  { 'on' 'off' }     'off';
                                'interp'      'string'  { 'on' 'off' }     'off';
                                'ersp'        'string'  { 'on' 'off' }     'off';
                                'recompute'   'string'  { 'on' 'off' }     'off';
                                'spec'        'string'  { 'on' 'off' }     'off';
                                'scalp'       'string'  { 'on' 'off' }     'off';
                                'allcomps'    'string'  { 'on' 'off' }     'off';
                                'itc'         'string'  { 'on' 'off' }     'off';
                                'savetrials'  'string'  { 'on','off' }     'off';
                                'rmicacomps'  'string'  { 'on' 'off' }     'off';
                                'design'      'integer' []                 STUDY.currentdesign;
                                'rmclust'     'integer' []                 [];
                                'rmbase'      'integer' []                 [];
                                'specparams'        'cell'    {}                 {};
                                'erspparams'        'cell'    {}                 {}}, 'std_precomp');
    if isstr(g), error(g); end;
    if ~isempty(g.rmbase), g.erpparams = { g.erpparams{:} 'rmbase' g.rmbase }; end;
    
    %%%%%JGwin modifications
    local_tw_flag = 0;
    local_twms_flag = 0;
    if isfield(g,'erspparams')
        if any(strcmp(g.erspparams,'subject tw matrix'))
            for i = find(strcmp(g.erspparams,'subject tw matrix'))
                if strcmp(g.erspparams{i-1},'timewarp')
                    local_tw_flag = 1;
                elseif strcmp(g.erspparams{i-1},'timewarpms')
                    local_twms_flag = 1;
                end
            end
        end
    end
    %%%
%commented out 12/7/2012    
% %     if isfield(g,'erspparams')
% %         if any(strcmp(g.erspparams,'subject tw matrix'))
% %             local_tw_flag = 1;
% %         end
% %     end
    %%%% 
    
    % union of all channel structures
    % -------------------------------
    computewhat = 'channels';
    if isstr(chanlist)
        if strcmpi(chanlist, 'channels')
            chanlist = [];
        else % components
            computewhat = 'components';
            if strcmpi(g.allcomps, 'on')
                chanlist = {};
                for index = 1:length(STUDY.datasetinfo)
                    chanlist = { chanlist{:} [1:size(ALLEEG(STUDY.datasetinfo(index).index).icaweights,1)] };
                end;
            else
                chanlist = { STUDY.datasetinfo.comps };
            end;
        end;
    end;
    if isempty(chanlist)
        alllocs = eeg_mergelocs(ALLEEG(:).chanlocs);
        chanlist = { alllocs.labels };
    elseif ~isnumeric(chanlist{1})
        alllocs = eeg_mergelocs(ALLEEG(:).chanlocs);
        [tmp c1 c2] = intersect( lower({ alllocs.labels }), lower(chanlist));
        [tmp c2] = sort(c2);
        alllocs = alllocs(c1(c2));
    end;
    
    % test if interp and reconstruct channel list
    % -------------------------------------------
    if strcmpi(computewhat, 'channels')
        if strcmpi(g.interp, 'on')
            STUDY.changrp = [];
            STUDY = std_changroup(STUDY, ALLEEG, chanlist, 'interp');
            g.interplocs = alllocs;
        else
            STUDY.changrp = [];
            STUDY = std_changroup(STUDY, ALLEEG, chanlist);
            g.interplocs = struct([]);
        end;
    end;
    
    % components or channels
    % ----------------------
    if strcmpi(computewhat, 'channels')
         curstruct = STUDY.changrp;
    else curstruct = STUDY.cluster;
    end;
    
    % compute ERPs
    % ------------
    if strcmpi(g.erp, 'on')
        % check dataset consistency
        % -------------------------
        allPnts = [ALLEEG([STUDY.design(g.design).cell.dataset]).pnts];
        if iscell(allPnts), allPnts = [ allPnts{:} ]; end;
        if length(unique(allPnts)) > 1
            error([ 'Cannot compute ERPs because datasets' 10 'do not have the same number of data points' ])
        end;
        
    for index = 1:length(STUDY.design(g.design).cell)
            desset = STUDY.design(g.design).cell(index);
            addopts = { 'savetrials', g.savetrials, 'recompute', g.recompute, 'fileout', desset.filebase, 'trialindices', desset.trials };
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_erp(ALLEEG(desset.dataset), 'channels', tmpchanlist, opts{:}, addopts{:}, g.erpparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_erp(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, addopts{:}, g.erpparams{:});
            end;
        end;
        if isfield(curstruct, 'erpdata')
            curstruct = rmfield(curstruct, 'erpdata');
            curstruct = rmfield(curstruct, 'erptimes');
        end;
    end;

    % compute component scalp maps
    % ----------------------------
    if strcmpi(g.scalp, 'on')
         %%%JGwin
        fprintf('Computing/checking topo files');
        %%%
        for index = 1:length(STUDY.datasetinfo)
            
            % find duplicate
            % --------------
            found = [];
            ind1 = STUDY.datasetinfo(index).index;
            inds = strmatch(STUDY.datasetinfo(index).subject, { STUDY.datasetinfo(1:index-1).subject });
            for index2 = inds'
                ind2 = STUDY.datasetinfo(index2).index;
                if isequal(ALLEEG(ind1).icawinv, ALLEEG(ind2).icawinv)
                    found = ind2;
                end;
            end;
            
            % make link if duplicate
            % ----------------------
            %%%JGwin
            %fprintf('Computing/checking topo file for dataset %d\n', ind1);
            fprintf('.');
            %%%
            if ~isempty(found)
                tmpfile1 = fullfile( ALLEEG(index).filepath, [ ALLEEG(index).filename(1:end-3) 'icatopo' ]); 
                tmp.file = fullfile( ALLEEG(found).filepath, [ ALLEEG(found).filename(1:end-3) 'icatopo' ]); 
                std_savedat(tmpfile1, tmp);
            else
                std_topo(ALLEEG(index), chanlist{index}, 'none', 'recompute', g.recompute);
            end;
        end;
        if isfield(curstruct, 'topo')
            curstruct = rmfield(curstruct, 'topo');
            curstruct = rmfield(curstruct, 'topox');
            curstruct = rmfield(curstruct, 'topoy');
            curstruct = rmfield(curstruct, 'topoall');
            curstruct = rmfield(curstruct, 'topopol');
        end;
         %%%JGwin
        fprintf(1,'%s\n','Done');
        %%%
    end;
    
    % compute spectrum
    % ----------------
    if strcmpi(g.spec, 'on')
        %%%JGwin - added this status update
        fprintf(1,'%s','Computing/checking spec files');
        %%%
        for index = 1:length(STUDY.design(g.design).cell)
             %%%JGwin - added status update
            fprintf('.');
            %%%
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_spec(ALLEEG(desset.dataset), 'channels', tmpchanlist, opts{:}, addopts{:}, g.specparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_spec(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, addopts{:}, g.specparams{:});
            end;
        end;
        if isfield(curstruct, 'specdata')
            curstruct = rmfield(curstruct, 'specdata');
            curstruct = rmfield(curstruct, 'specfreqs');
        end;
        %%%JGwin - added status update
        fprintf(1,'%s\n','Done');
        %%%
        
    end;

    % compute ERSP and ITC
    % --------------------
    if strcmpi(g.ersp, 'on') | strcmpi(g.itc, 'on')
        % check dataset consistency
        % -------------------------
        allPnts = [ALLEEG([STUDY.design(g.design).cell.dataset]).pnts];
        if iscell(allPnts), allPnts = [ allPnts{:} ]; end;
        if length(unique(allPnts)) > 1
            error([ 'Cannot compute ERSPs/ITCs because datasets' 10 'do not have the same number of data points' ])
        end;
        
        if strcmpi(g.ersp, 'on') & strcmpi(g.itc, 'on'), type = 'both';
        elseif strcmpi(g.ersp, 'on')                   , type = 'ersp';
        else                                             type = 'itc';
        end;
        
        % check for existing files
        % ------------------------
        %%%%JGwin
        if strcmpi(g.recompute, 'off')
            guimode = 'usedisk';
        else
            guimode = 'recompute';
        end
        %guimode = 'guion';
        %%%%
        [ tmpX tmpt tmpf g.erspparams ] = std_ersp(ALLEEG(1), 'channels', 1, 'type', type, 'recompute', 'on', 'getparams', 'on', 'savetrials', g.savetrials, g.erspparams{:});
        if strcmpi(g.recompute, 'off')
            for index = 1:length(STUDY.design(g.design).cell)
                desset = STUDY.design(g.design).cell(index);
                if strcmpi(computewhat, 'channels')
                     filename = [ desset.filebase '.datersp'];
                else filename = [ desset.filebase '.icaersp'];
                end;
                %%%JGwin: ignore timewarp params if using subject specific
                %%%timewarp values
                if any(strcmp(g.erspparams,'subject tw matrix'))
                    ignore_fields = { 'plotitc' 'plotersp' 'plotphase' 'savetrials' 'timewarp' 'timewarpms'};
                    %don't want to change the timewarp fields so don't get
                    %new ersp params
                    [guimode ] = std_filecheck(filename, g.erspparams, guimode, ignore_fields);
                else
                    ignore_fields = { 'plotitc' 'plotersp' 'plotphase' 'savetrials'};
                    [guimode, g.erspparams] = std_filecheck(filename, g.erspparams, guimode, ignore_fields);
                end
                %%%
                if strcmpi(guimode, 'cancel'), return; end;

             end;
            if strcmpi(guimode, 'usedisk') | strcmpi(guimode, 'same'), g.recompute = 'off'; 
            else                                                       g.recompute = 'on'; 
            end;
        end;
        
        % check for existing files
        % ------------------------
        if isempty(g.erspparams), 
            tmpparams = {};
        elseif iscell(g.erspparams), 
            tmpparams = g.erspparams; 
        else
            tmpparams      = fieldnames(g.erspparams); tmpparams = tmpparams';
            tmpparams(2,:) = struct2cell(g.erspparams);
        end;
        tmpparams = { tmpparams{:} 'recompute' g.recompute };
        
        for index = 1:length(STUDY.design(g.design).cell) %change 3 back to 1
            desset = STUDY.design(g.design).cell(index);
            %%%%JGwin modification to use tw saved in ALLEEG
            tw_idx = find(strcmp(tmpparams,'timewarp'));
            if ~isempty(tw_idx) && local_tw_flag
                tmpparams(tw_idx+1) = {ALLEEG(desset.dataset).timewarp.latencies};
            end
            twms_idx = find(strcmp(tmpparams,'timewarpms'));
            if ~isempty(twms_idx) && local_twms_flag
                tmpparams(twms_idx+1) = {ALLEEG(desset.dataset).timewarp.latencies};
            end
            %%%%
            
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_ersp(ALLEEG(desset.dataset), 'channels', tmpchanlist, 'type', type, 'fileout', desset.filebase, 'trialindices', desset.trials, opts{:}, tmpparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_ersp(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, 'type', type, 'fileout', desset.filebase, 'trialindices', desset.trials, tmpparams{:});
            end;
        end;
        if isfield(curstruct, 'erspdata')
            curstruct = rmfield(curstruct, 'erspdata');
            curstruct = rmfield(curstruct, 'ersptimes');
            curstruct = rmfield(curstruct, 'erspfreqs');
        end;
        if isfield(curstruct, 'itcdata')
            curstruct = rmfield(curstruct, 'itcdata');
            curstruct = rmfield(curstruct, 'itctimes');
            curstruct = rmfield(curstruct, 'itcfreqs');
        end;
    end;

    % components or channels
    % ----------------------
    if strcmpi(computewhat, 'channels')
         STUDY.changrp = curstruct;
    else STUDY.cluster = curstruct;
    end;
    
    return;

%     % get channel indices from changrp structure
%     % ------------------------------------------    
%     function chaninds = getchannelindices(changrp, datasetind);    
%     
%     chaninds = [];
%     for index = 1:length(changrp)    
%         tmpind = find([ changrp(index).setinds{:} ] == datasetind);
%         if ~isempty(tmpind)
%             chaninds = [ chaninds index ];
%         end;
%     end;
        
    % find components in cluster for specific dataset
    % -----------------------------------------------
    function rmcomps = getclustcomps(STUDY, rmclust, settmpind);    
    
        rmcomps   = cell(1,length(settmpind));
        for idat = 1:length(settmpind) % scan dataset for which to find component clusters
            currentdat = settmpind(idat);
            for rmi = 1:length(rmclust) % scan clusters
                compinds = STUDY.cluster(rmclust(rmi)).allinds;
                
                for ind = 1:length(STUDY.cluster(rmclust(rmi)).setinds(:)) % scan datasets of cluster
                    setinds  = [STUDY.design(STUDY.currentdesign).cell(STUDY.cluster(rmclust(rmi)).setinds{ind}).dataset ];
                    indmatch = find(setinds == currentdat);
                    if ~isempty(indmatch)
                        rmcomps{idat} = union( rmcomps{idat}, compinds{ind}(indmatch));
                    end;
                end;
            end;
        end;
        
    % make option array and channel list (which depend on interp) for any type of measure
    % ----------------------------------------------------------------------
    function [tmpchanlist, opts] = getchansandopts(STUDY, ALLEEG, chanlist, idat, g);
        
        opts = { };
        if ~isempty(g.rmclust)
            opts = { opts{:} 'rmcomps' getclustcomps(STUDY, g.rmclust, idat) };                
        elseif strcmpi(g.rmicacomps, 'on')
            for ind = 1:length(idat)
                rmcomps{ind} = find(ALLEEG(idat(1)).reject.gcompreject);
            end;
            opts = { opts{:} 'rmcomps' rmcomps };
        end;
        if strcmpi(g.interp, 'on')
            tmpchanlist = chanlist;
            allocs = eeg_mergelocs(ALLEEG.chanlocs);
            [tmp1 tmp2 neworder] = intersect( {allocs.labels}, chanlist);
            [tmp1 ordertmp2] = sort(tmp2);
            neworder = neworder(ordertmp2);
            opts = { opts{:} 'interp' allocs(neworder) };
        else
            newchanlist = [];
            tmpchanlocs = ALLEEG(idat(1)).chanlocs;
            chanlocs = { tmpchanlocs.labels };
            for i=1:length(chanlist)
                newchanlist = [ newchanlist strmatch(chanlist{i}, chanlocs, 'exact') ];
            end;
            tmpchanlocs =  ALLEEG(idat(1)).chanlocs;
            tmpchanlist = { tmpchanlocs(newchanlist).labels };
        end;