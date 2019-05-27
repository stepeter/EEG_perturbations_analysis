% mod_pop_runemd() - Run an EMD decomposition of an EEG dataset using runemd(), 
%                runemdtrialbytrial(), memd(), weighted_slidingemd() 
%               
% Usage:
%   >> OUT_EEG = mod_pop_runemd( EEG ); % pops-up a data entry window
%   >> OUT_EEG = mod_pop_runemd( EEG, 'key', 'val' ); % no pop_up
%
% Graphic interface:
%   "EMD algorithm to use" - [pop-up menue] The EMD algorithm to use for 
%                 EMD decomposition. Command line equivalent: 'emdtype'
%   "Commandline options" - [edit box] Command line options to forward
%                 to the EMD algorithm. Command line equivalent: 'options' 
% Inputs:
%   EEG         - input EEG dataset or array of datasets
%
% Optional inputs:
%   'emdtype'    - ['EMD'|'EEMD'|'MEMD'|'wSEMD']  algorithm 
%                 to use for the EMD decomposition. The nature of any 
%                 differences in the results of these algorithms have 
%                 not been well characterized.
%   'norm'        - ['on'|'off'] 'on' normalize input dataset 
%   'trialbytrial' - ['on'|'off'] 'on' decompose input dataset trial by trial
%   'dataset'   - [integer array] dataset index or indices.
%   'nmodes'    - [integer number] number of modes (in case of EMD and EEMD).
%   'ensemblenum'    - [integer number] number of ensemble (in case of EEMD).
%   'noiseassist'    - [integer number] number of assisted noise (in case of EEMD).
%   'wsize'    - [integer number] window size (in case of wSEMD).
%   'ssize'    - [integer number] step size (in case of wSEMD).
%   'chanind'   - [integer array or cell array] subset of channel indices 
%                 for running the EMD decomposition. Alternatively, you may
%                 also enter channel types here in a cell array.
%   
%   'key','val' - EMD algorithm options (see EMD routine help messages).
% 
%
% Outputs:
%   OUT_EEG = The input EEGLAB dataset with new fields IMFs, imfnumber, eemdensemblenumber,
%   eemdassistednoise, emdchansind, emdchanlocs, emdnbchan.
%
%    See also: runemd(),runemdtrialbytrial(),memd() and weighted_slidingemd().
%
%
%
% Acknowledgment: This plugin is basically based on the EEGLAB Toolbox code, publicly available from
%                http://sccn.ucsd.edu/eeglab/downloadtoolbox.html , EMD code prepared by Zhaohua Wu (zhwu@cola.iges.org)
%                ,MEMD prepared by Naveed ur Rehman and Danilo P. Mandic, Oct-2009 and wSEMD prepared by  Angela Zeiler  
%                (angela.zeiler@biologie.uni-regensburg.de). 
%
%
% Authors:Karema Al-Subari (Karema.Al-Subari@ur.de) and Saad Al-Baddai (Saad.Al-Baddai@ur.de), CIML group,Regensburg Uni, 2014


% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License of EEGLAB as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%




function [IMFs, com] = mod_pop_runemd( EEG, varargin )

com = '';

 


if nargin < 1   
    help mod_pop_runemd;
    return;
end;

if isempty(EEG.data)
     errordlg2(strvcat('Error: no data to decompose. use menu "File > Import or Load Data " first.'), 'Error');
          return;
 end;


% find available algorithms
% -------------------------
allalgs   = { 'EMD' 'EEMD' 'MEMD' 'wSEMD'}; % do not use egld_EMD => too slow
selectalg = {};
linenb    = 1;
count     = 1;


% special AMEMD
% -------------
selectamemd = 0;
winsize=100;%%windowsize
stsize=20;%%stepsize
defaultensemble = [ '2' ] ;
defaultnoise = [ '0.1' ] ;
defaultmodenumber = [ '4' ] ;
defaultnoisechannel = 4;
if nargin > 1
    if isstr(varargin{1})
        if strcmpi(varargin{1}, 'selectamemd')
            selectamemd = 1;
            allalgs = { 'amemd' allalgs{:} };
            defaultopts = sprintf('''outdir'', ''%s''', fullfile(pwd, 'amemdout'));
        elseif strcmpi(varargin{1}, 'selectamemdloc')
            selectamemd = 1;
            allalgs = { 'amemd' allalgs{:} };
            defaultopts = sprintf('''outdir'', ''%s'', ''qsub'', ''off''', fullfile(pwd, 'amemdout'));
        end;
    end;
end;



% popup window parameters
% -----------------------
fig = [];

if nargin < 2 | selectamemd
    commandchans = [ 'tmpchans = get(gcbf, ''userdata'');' ...
                     'tmpchans = tmpchans{1};' ...
                     'set(findobj(gcbf, ''tag'', ''chantype''), ''string'', ' ...
                     '       int2str(pop_chansel( tmpchans )));' ...
                     'clear tmpchans;' ];

   
    cb_emd = [ 'if get(gcbo, ''value'') > 3, ' ...
               '     set(findobj(gcbf, ''tag'', ''params''), ''string'', ''''''extended'''', 1'');' ...
               'else set(findobj(gcbf, ''tag'', ''params''), ''string'', '''');' ...
               'end;'.... % ui handle
               'E= findobj(''tag'',''params4'',''parent'',gcbf);'....
               'items = get(E, ''value'');'...
               ' if items==1 & EEG.trials==1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''off'');' ...
               ' elseif items==1 & EEG.trials~=1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''off'');' ...
               ' elseif items==2 & EEG.trials==1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''off'');' ...
               ' elseif items==2 & EEG.trials~=1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''off'');' ...
               ' elseif items==3 & EEG.trials==1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''on'');' ...
               ' elseif items==3 & EEG.trials~=1;' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''on'');' ...
               'else,' ...
               ' set(findobj(gcf, ''tag'', ''params1''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params2''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params3''), ''enable'', ''off'');' ...
               ' set(findobj(gcf, ''tag'', ''params6''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params7''), ''enable'', ''on'');' ...
               ' set(findobj(gcf, ''tag'', ''params8''), ''enable'', ''off'');' ...
               ' end;'
               
                 ];
 
  if EEG.trials==1
    promptstr    = { { 'style' 'text'       'string' 'EMD algorithm to use (click to select)' 'fontweight', 'bold'} ...
                     { 'style' 'popup'    'string' strvcat(allalgs{:}) 'tag' 'params4' 'callback', cb_emd } ...
                     { 'style' 'checkbox'       'string' 'Normalize' } ...
                     { 'style' 'text'       'string' 'single-trial Decomposition' 'tag' 'params5'  } ...
                     { 'style' 'text'       'string' 'number of modes ' } ...
                     { 'style' 'text'       'string' 'ensemble number ' } ...
                     { 'style' 'text'       'string' 'noise-assisted ' } ...                    
                     { 'style' 'edit'       'string' defaultmodenumber 'tag' 'params1' } ...
                     { 'style' 'edit'       'string' defaultensemble 'tag' 'params2' 'enable' 'off' } ...
                     { 'style' 'edit'       'string' defaultnoise 'tag' 'params3' 'enable' 'off' } ...
                     { 'style' 'text'       'string' 'window size' } ...
                     { 'style' 'text'       'string' 'step size ' } ... 
                     { 'style' 'text'       'string' 'number of noise channels ' } ... 
                     { 'style' 'edit'       'string' winsize 'tag' 'params6' 'enable' 'off' } ...
                     { 'style' 'edit'       'string' stsize 'tag' 'params7' 'enable' 'off' } ...
                     { 'style' 'edit'       'string' defaultnoisechannel 'tag' 'params8' 'enable' 'off' } ...
                     { 'style' 'text'       'string' 'Channel type(s) or channel indices' } ...
                     { 'style' 'edit'       'string' '' 'tag' 'chantype' }  ...
                     { 'style' 'pushbutton' 'string' '... channels' 'callback' commandchans } };
  else
  promptstr    = { { 'style' 'text'       'string' 'EMD algorithm to use (click to select)' 'fontweight', 'bold' } ...
                     { 'style' 'popup'    'string' strvcat(allalgs{:}) 'tag' 'params4' 'callback', cb_emd } ...
                     { 'style' 'checkbox'       'string' 'Normalize' } ...
                     { 'style' 'text'       'string' 'Trial-by-Trial Decomposition ' 'tag' 'params5' } ...
                     { 'style' 'text'       'string' 'number of modes ' } ...
                     { 'style' 'text'       'string' 'ensemble number ' } ...
                     { 'style' 'text'       'string' 'noise-assisted ' } ...                    
                     { 'style' 'edit'       'string' defaultmodenumber 'tag' 'params1' } ...
                     { 'style' 'edit'       'string' defaultensemble 'tag' 'params2' 'enable' 'off'} ...                    
                     { 'style' 'edit'       'string' defaultnoise 'tag' 'params3' 'enable' 'off' } ...    
                     { 'style' 'text'       'string' 'window size' } ...
                     { 'style' 'text'       'string' 'step size ' } ...  
                     { 'style' 'text'       'string' 'number of noise channels ' } ...                 
                     { 'style' 'edit'       'string' winsize 'tag' 'params6' 'enable' 'off' } ...                    
                     { 'style' 'edit'       'string' stsize 'tag' 'params7' 'enable' 'off' } ...
                     { 'style' 'edit'       'string' defaultnoisechannel 'tag' 'params8' 'enable' 'off' } ...
                     { 'style' 'text'       'string' 'Channel type(s) or channel indices' } ...
                     { 'style' 'edit'       'string' '' 'tag' 'chantype' }  ...
                     { 'style' 'pushbutton' 'string' '... channels' 'callback' commandchans } };
  end

    geometry = { [0.8 0.4 0.4 0.5] [] [0.5 0.5 0.5] [0.5 0.5 0.5] [] [0.5 0.5 0.5] [0.5 0.5 0.5] [] [1 1 0.5] };

  
                     
    % channel types
    % -------------
    if isfield(EEG.chanlocs, 'type'), 
        tmpchanlocs = EEG.chanlocs;
        alltypes = { tmpchanlocs.type };
        indempty = cellfun('isempty', alltypes);
        alltypes(indempty) = '';
        try, 
            alltypes = unique(alltypes);
        catch, 
            alltypes = '';
        end;
    else
        alltypes = '';
    end;
    
    % channel labels
    % --------------
    if ~isempty(EEG.chanlocs)
        tmpchanlocs = EEG.chanlocs;        
           alllabels = { tmpchanlocs.labels };
           
    else
        for index = 1:EEG.nbchan
            alllabels{index} = int2str(index);
        end;
    end;
    
    % gui
    % ---
    result       = inputgui( 'geometry', geometry, 'uilist', promptstr, ...
                             'helpcom', 'pophelp(''mod_pop_runemd'')', ...
                             'title', 'Run EMD decomposition -- mod_pop_runemd()', 'userdata', { alllabels alltypes } );

    if length(result) == 0 return; end;        
    options = { 'emdtype' allalgs{result{1}}  'norm' (result{2}) 'nmodes' eval(result{3}) 'ensemblenum' eval(result{4}) 'noiseassist' eval(result{5}) 'wsize' eval(result{6}) 'ssize' eval(result{7}) 'defaultnoisechannel' eval(result{8}) 'dataset' [1:length(EEG)]  };
    if ~isempty(result{9})
        if ~isempty(str2num(result{9})), options = { options{:} 'chanind' str2num(result{9}) };
        else                             options = { options{:} 'chanind' parsetxt(result{9}) }; 
        end;
    end;

   
else 
    if mod(length(varargin),2) == 1
        options = { 'emdtype' varargin{1:end} };
    else
        options = varargin;
    end;
end;

% decode input arguments
% ----------------------
[ g addoptions ] = finputcheck( options, { 'emdtype'        'string'  allalgs   'runemd'; ...
                            'dataset'        'integer' []        [1:length(EEG)];
                            'options'        'cell'    []        {};
                            'concatenate'    'string'  { 'on','off' }   'off';
                            'concatcond'     'string'  { 'on','off' }   'off';
                            'chanind'        { 'cell','integer' } { [] [] }        [];}, ...
                            'mod_pop_runemd', 'ignore');
if isstr(g), error(g); end;
if ~isempty(addoptions), g.options = { g.options{:} addoptions{:}}; end;


% select datasets, create new big dataset if necessary
% ----------------------------------------------------
if length(g.dataset) == 1
    EEG = EEG(g.dataset);
elseif length(EEG) > 1 & ~strcmpi(g.concatenate, 'on') & ~strcmpi(g.concatcond, 'on')
    [ EEG com ] = eeg_eval( 'mod_pop_runemd', EEG, 'warning', 'off', 'params', ...
           { 'emdtype' g.emdtype 'options' g.options 'chanind' g.chanind } );
    return;
elseif length(EEG) > 1 & strcmpi(g.concatcond, 'on')
    allsubjects = { EEG.subject };
    allsessions = { EEG.session };
    allgroups   = { EEG.group };
    alltags     = zeros(1,length(allsubjects));
    if any(cellfun('isempty', allsubjects))
        disp('Aborting: Subject names missing from at least one dataset.');
        return;
    end;
    dats = {};
    for index = 1:length(allsubjects)
        if ~alltags(index)
            allinds = strmatch(allsubjects{index}, allsubjects, 'exact');
            rmind = [];
            for tmpi = 2:length(allinds)
                if ~isequal(allsessions(allinds(1)), allsessions(allinds(tmpi))), rmind = [rmind tmpi];
                elseif ~isequal(allgroups(allinds(1)), allgroups(allinds(tmpi))), rmind = [rmind tmpi]; 
                end;
            end;
            allinds(rmind) = [];
            fprintf('Found %d datasets for subject ''%s''\n', length(allinds), allsubjects{index});
            dats = { dats{:} allinds };
            alltags(allinds) = 1;
        end;
    end;
 
    fprintf('**************************\nNOW RUNNING ALL DECOMPOSITIONS\n****************************\n');
    for index = 1:length(dats)
        EEG(dats{index}) = mod_pop_runemd(EEG(dats{index}), 'emdtype', g.emdtype, ...
            'options', g.options, 'chanind', g.chanind, 'concatenate', 'on');
        for idat = 1:length(dats{index})
            EEG(dats{index}(idat)).saved = 'off';
            pop_saveset(EEG(dats{index}(idat)), 'savemode', 'resave');
            EEG(dats{index}(idat)).saved = 'on';
        end;
    end;
    com = sprintf('%s = mod_pop_runemd(%s, %s);', inputname(1),inputname(1), ...
              vararg2str({ 'emdtype' g.emdtype 'concatcond' 'on' 'options' g.options }) );
    return;
else
    disp('Concatenating datasets...');
    EEG = EEG(g.dataset(1));
    
    % compute total data size
    % -----------------------
    totalpnts = 0;
    for i = g.dataset
        totalpnts = totalpnts+EEG(g.dataset(i)).pnts*EEG(g.dataset(i)).trials;
    end;
    EEG.data = zeros(EEG.nbchan, totalpnts);
    
    % copy data
    % ---------
    cpnts = 1;
    for i = g.dataset
        tmplen = EEG(g.dataset(i)).pnts*EEG(g.dataset(i)).trials;
        TMP = eeg_checksetemd(EEG(g.dataset(i)), 'loaddata');
        EEG.data(:,cpnts:cpnts+tmplen-1) = reshape(TMP.data, size(TMP.data,1), size(TMP.data,2)*size(TMP.data,3));
        cpnts = cpnts+tmplen;
    end;
   
    EEG.trials = 1;
    EEG.pnts   = size(EEG.data,2);
end;    

% Store and then remove current EEG emd weights and sphere

    fprintf('Saving current IMFs decomposition in "EEG.IMFs" (etc.).\n');
   
    


%% select all chanels

  if isempty(g.chanind)
     g.chanind =(1:EEG.nbchan);
  end

%   if isstruct(g.chanind) 
% 
%     g.chanid = (1:EEG.nbchan);
if iscell(g.chanind)
    g.chanind = eeg_chantype(EEG.chanlocs, g.chanind);
end;




%------------------------------




%------------------------------
% Normalize dataset
% -----------------------------
if EEG.trials==1
 if g.norm==0
  tmpdata=reshape( EEG.data(g.chanind,:,:), length(g.chanind), EEG.pnts*EEG.trials);
 else
  tmpdata = reshape( EEG.data(g.chanind,:,:), length(g.chanind), EEG.pnts*EEG.trials);
  %% normalize data zscore

  tmpdata = (zscore_normalize(tmpdata'))'; % zero mean 
 end
else
 if g.norm==0
  tmpdata=EEG.data(g.chanind,:,:);
 else
   tmpdata =  EEG.data(g.chanind,:,:);
  
   %%tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean
  
   %% normalize data zscore
   for ii = 1: length(g.chanind)
     tmpdataNormalize(ii,:,:) = zscore_normalize((tmpdata(ii,:,:)));
   end
   tmpdata = tmpdataNormalize;
   
  end

end
  
switch lower(g.emdtype)
  
    case 'emd' 
        try, if ismatlab, g.options = {  g.options{:}, 'interupt', 'on' }; end; catch, end; 
          if ~isempty(tmpdata)

             if EEG.trials==1
                IMFs = runemd(tmpdata,'Nstd',0,'NE',1,'modes',g.nmodes,'trials',EEG.trials);
             else
                IMFs = runemdtrialbytrial(tmpdata,'Nstd',0,'NE',1,'modes',g.nmodes,'trials',EEG.trials);
             end
          EEG.imfnumber=g.nmodes+1;
          EEG.emdchansind=g.chanind;
          EEG.emdchanlocs=EEG.chanlocs(g.chanind);
          EEG.emdnbchan=length(g.chanind); 
                

       
         else 
            disp(['There is no data to decompose!']);
           
         end;
     case 'eemd'
         if ~isempty(tmpdata)
           if EEG.trials==1
             
             IMFs = runemd(tmpdata,'Nstd',g.noiseassist,'NE',g.ensemblenum,'modes',g.nmodes,'trials',EEG.trials); 
            else
             
             IMFs = runemdtrialbytrial(tmpdata,'Nstd',g.noiseassist,'NE',g.ensemblenum,'modes',g.nmodes,'trials',EEG.trials); 
            end
           EEG.imfnumber=g.nmodes+1;
           EEG.eemdensemblenumber=g.ensemblenum;
           EEG.eemdassistednoise=g.noiseassist;
           EEG.emdchansind=g.chanind;  
           EEG.emdchanlocs=EEG.chanlocs(g.chanind);
           EEG.emdnbchan=length(g.chanind); 
     

        else 
            disp(['There is no data to decompose!']);
           
        end;

   case 'memd'
     tmpdat.data=double(tmpdata);
%     tmpdat.check=g.trialbytrial;
     tmpdat.trials=EEG.trials;
 
         if ~isempty(tmpdat)
            if EEG.trials==1
             IMFs = ena_memd(tmpdat,g.ensemblenum,g.defaultnoisechannel,g.noiseassist); 
            else
             IMFs = ena_memdtrialbytrial(tmpdat,g.ensemblenum,g.defaultnoisechannel,g.noiseassist); 
  
            end
            
                
             EEG.imfnumber=size(IMFs,2);
             EEG.emdchansind=g.chanind;  
             EEG.emdchanlocs=EEG.chanlocs(g.chanind);
             EEG.emdnbchan=length(g.chanind); 
          else 
            disp(['There is no data to decompose!']);
           
        end;
        
        
        case 'wsemd' 
        try, if ismatlab, g.options = {  g.options{:}, 'interupt', 'on' }; end; catch, end; 
          if ~isempty(tmpdata)
               if EEG.trials==1
                   
                IMFs = weighted_slidingemd(tmpdata','windowsize',g.wsize,'stepsize',g.ssize,'ntrials',EEG.trials);
               else
                   
                 error('Error: SEMD works only in case there is only one epoche!');
                 return;
               end
                EEG.imfnumber=g.nmodes+1;
                EEG.emdchansind=g.chanind;
                EEG.emdchanlocs=EEG.chanlocs(g.chanind);
                EEG.emdnbchan=length(g.chanind); 
                

       
         else 
            disp(['There is no data to decompose!']);
           
         end;
        
        
    

        clear tmp;
        close(fig);
     otherwise, error('mod_pop_runemd: unrecognized algorithm');
end;
if ~isempty(fig), try, close(fig); catch, end; end;
% copy back data to datasets if necessary
% ---------------------------------------
if length(g.dataset) > 1
            
    EEG = eeg_checksetemd(EEG);
else

    EEG = eeg_checksetemd(EEG);
    EEG = eeg_store(EEG, EEG, g.dataset);

end
if nargin < 2 || selectamemd
    com = sprintf('%s = mod_pop_runemd(%s, %s);', inputname(1), inputname(1),  vararg2str(options) ); %vararg2str({ 'emdtype' g.emdtype 'dataset' g.dataset 'options' g.options }) );
end;



W_MAIN = findobj('tag', 'EEGLAB');
%file_m = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'File');  set(file_m, 'enable', 'on');
EMDmenu = findobj('parent', W_MAIN, 'type', 'uimenu', 'label', 'EMDLAB');
temd_m = findobj('parent', EMDmenu, 'type', 'uimenu', 'label', 'ERM Maps');
ERMC_m = findobj('parent', EMDmenu, 'type', 'uimenu', 'label', 'ERM and Maps');




    if EEG.trials == 1 & isfield(EEG,'IMFs')
			

            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Scrolling Modes')   , 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Power Spectra and Maps')        , 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Mode Properties'), 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
  end    


 % epoched data

    if EEG.trials>1 & isfield(EEG,'IMFs')
  	
            set( findobj('parent', EMDmenu, 'type', 'uimenu','Label','Scrolling Modes'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Power Spectra and Maps'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'ERM Maps'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', temd_m, 'type', 'uimenu', 'Label', 'In 2-D'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', temd_m, 'type', 'uimenu', 'Label', 'In 3-D'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Mode Properties'), 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Hilbert-Huang/Fourier Transform')   , 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'IMFs and ERM')   , 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'ERM and Maps')   , 'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on'); 
            set( findobj('parent', ERMC_m, 'type', 'uimenu', 'Label', 'With modes maps'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', ERMC_m, 'type', 'uimenu', 'Label', 'In rectangular array'),'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
            set( findobj('parent', EMDmenu, 'type', 'uimenu', 'Label', 'Compare ERMs')   ,'userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');% 
    end
      
return;
          
            
