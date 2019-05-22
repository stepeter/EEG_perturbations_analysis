function [selection,value] = eegchan_listdlg(EEG,varargin)
%EEGCHAN_LISTDLG  List selection dialog box.
%   [SELECTION,OK] = EEGCHAN_LISTDLG(EEG,'ListString',S) creates a modal dialog box
%   which allows you to select a string or multiple channels from a list.
%   An eegplot and topoplot are created from the input var EEG. Channels
%   selected in the listbox are highlighted in the corresponding eegplot and
%   topoplot. This function is based on Matlab's listdlg box.
%   
%   eegchan_listdlg created by HJH, 04-19-13   
%
%   SELECTION is a vector of indices of the selected strings (length 1 in
%   the single selection mode).  This will be [] when OK is 0.  OK is 1 if
%   you push the OK button, or 0 if you push the Cancel button or close the
%   figure.
%
%   Double-clicking on an item or pressing <CR> when multiple items are
%   selected has the same effect as clicking the OK button.  Pressing <CR>
%   is the same as clicking the OK button. Pressing <ESC> is the same as
%   clicking the Cancel button.
%   
%
%   Parameter inputs are in parameter,value pairs:
%
%   Parameter       Description
%   'ListString'    cell array of strings for the list box.
%   'SelectionMode' string; can be 'single' or 'multiple'; defaults to
%                   'multiple'.
%   'ListSize'      [width height] of listbox in pixels; defaults
%                   to [160 300].
%   'InitialValue'  vector of indices of which items of the list box
%                   are initially selected; defaults to the first item.
%   'Name'          String for the figure's title; defaults to ''.
%   'PromptString'  string matrix or cell array of strings which appears 
%                   as text above the list box; defaults to {}.
%   'OKString'      string for the OK button; defaults to 'OK'.
%   'CancelString'  string for the Cancel button; defaults to 'Cancel'.
%
%   A 'Select all' button is provided in the multiple selection case.
%
%   Example:
%     d = dir;
%     str = {d.name};
%     [s,v] = listdlg('PromptString','Select a file:',...
%                     'SelectionMode','single',...
%                     'ListString',str)
 %
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG,
%    MSGBOX, QUESTDLG, WARNDLG.

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.20.4.11 $  $Date: 2010/04/21 21:34:03 $

%   'uh'            uicontrol button height, in pixels; default = 22.
%   'fus'           frame/uicontrol spacing, in pixels; default = 8.
%   'ffs'           frame/figure spacing, in pixels; default = 8.

% simple test:
%
% d = dir; [s,v] = listdlg('PromptString','Select a file:','ListString',{d.name});
% 

% Generate a warning in -nodisplay and -noFigureWindows mode.
% warnfiguredialog('listdlg');

error(nargchk(1,inf,nargin))

figname = '';
smode = 2;   % (multiple)
promptstring = {};
liststring = [];
listsize = [160 300];
initialvalue = [];
okstring = 'OK';
cancelstring = 'Cancel';
fus = 8;
ffs = 8;
uh = 22;

% create eegplot
eegplot(EEG.data,'command',[],'eloc_file',EEG.chanlocs,'winlength',50);
h1 = gcf;
h_eeg = findobj(h1,'tag','eegaxis');
h_eegchan = get(h_eeg,'children');

% create chan loc plot
h2 = figure;
mod_topoplot([],EEG.chanlocs,'electrodes','on','emarker',{[1:EEG.nbchan],'.','k',10,1});

h_topo = get(h2,'children');
hs = get(h_topo,'children');
h_elocs = hs(1);
x = get(h_elocs,'XData');
y = get(h_elocs,'YData');

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error('MATLAB:listdlg:InvalidArgument', 'Arguments to LISTDLG must come param/value in pairs.')
end
for i=1:2:length(varargin)
    switch lower(varargin{i})
     case 'name'
      figname = varargin{i+1};
     case 'promptstring'
      promptstring = varargin{i+1};
     case 'selectionmode'
      switch lower(varargin{i+1})
       case 'single'
        smode = 1;
       case 'multiple'
        smode = 2;
      end
     case 'listsize'
      listsize = varargin{i+1};
     case 'liststring'
      liststring = varargin{i+1};
     case 'initialvalue'
      initialvalue = varargin{i+1};
     case 'uh'
      uh = varargin{i+1};
     case 'fus'
      fus = varargin{i+1};
     case 'ffs'
      ffs = varargin{i+1};
     case 'okstring'
      okstring = varargin{i+1};
     case 'cancelstring'
      cancelstring = varargin{i+1};
     otherwise
      error('MATLAB:listdlg:UnknownParameter', 'Unknown parameter name passed to LISTDLG.  Name was %s', varargin{i})
    end
end

if ischar(promptstring)
    promptstring = cellstr(promptstring); 
end

if isempty(initialvalue)
    initialvalue = 1;
end

if isempty(liststring)
    error('MATLAB:listdlg:NeedParameter', 'ListString parameter is required.')
end

ex = get(0,'DefaultUicontrolFontSize')*1.7;  % height extent per line of uicontrol text (approx)

fp = get(0,'DefaultFigurePosition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'normal' ...
    'visible'                'off' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
            };

liststring=cellstr(liststring);

% initial channel selection status to be off, 0
chan_sel_status = zeros(1,length(liststring));
[sel_chan_names ia eeg_chan_num] = intersect(liststring,{EEG.chanlocs.labels});

fig = figure(fig_props{:});

if length(promptstring)>0
    prompt_text = uicontrol('style','text','string',promptstring,...
        'horizontalalignment','left',...
        'position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring)) ...
        listsize(1) ex*length(promptstring)]); %#ok
end

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;

listbox = uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
                    'string',liststring,...
                    'backgroundcolor','w',...
                    'max',smode,...
                    'tag','listbox',...
                    'value',initialvalue, ...
                    'callback', {@doListboxClick});

ok_btn = uicontrol('style','pushbutton',...
                   'string',okstring,...
                   'position',[ffs+fus ffs+fus btn_wid uh],...
                   'Tag','ok_btn',...
                   'callback',{@doOK,listbox, h1,h2});

cancel_btn = uicontrol('style','pushbutton',...
                       'string',cancelstring,...
                       'position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
                       'Tag','cancel_btn',...
                       'callback',{@doCancel,listbox, h1,h2});

if smode == 2
    selectall_btn = uicontrol('style','pushbutton',...
                              'string','Select all',...
                              'position',[ffs+fus 4*fus+ffs+uh listsize(1) uh],...
                              'tag','selectall_btn',...
                              'callback',{@doSelectAll, listbox, h_eeg, eeg_chan_num, h_topo,x,y});

    if length(initialvalue) == length(liststring)
        set(selectall_btn,'enable','off')
    end
    set(listbox,'callback',{@doListboxClick, selectall_btn, fig, h_eeg, chan_sel_status, eeg_chan_num, sel_chan_names, liststring, h_topo,x,y})
end

set([fig, ok_btn, cancel_btn, listbox], 'keypressfcn', {@doKeypress, listbox});

% set(fig,'position',getnicedialoglocation(fp, get(fig,'Units')));
% Make ok_btn the default button.
% setdefaultbutton(fig, ok_btn);

% make sure we are on screen
movegui(fig)
set(fig, 'visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(listbox);
    uiwait(fig);
catch
    if ishghandle(fig)
        delete(fig)
    end
end

if isappdata(0,'ListDialogAppData__')
    ad = getappdata(0,'ListDialogAppData__');
    selection = ad.selection;
    value = ad.value;
    rmappdata(0,'ListDialogAppData__')
else
    % figure was deleted
    selection = [];
    value = 0;
    delete(h1); delete(h2);
end


%% figure, OK and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
 case 'escape'
  doCancel([],[],listbox);
end

%% OK callback
function doOK(ok_btn, evd, listbox,h1,h2) %#ok
if (~isappdata(0, 'ListDialogAppData__'))
    ad.value = 1;
    ad.selection = get(listbox,'value');
    
    setappdata(0,'ListDialogAppData__',ad);
    delete(gcbf);
    delete(h1); delete(h2);
end

%% Cancel callback
function doCancel(cancel_btn, evd, listbox, h1,h2) %#ok
ad.value = 0;
ad.selection = [];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);
delete(h1); delete(h2);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox, h_eeg, eeg_chan_num, h_topo,x,y) %#ok
set(selectall_btn,'enable','off')
set(listbox,'value',1:length(get(listbox,'string')));
selected_chan = get(listbox,'value');
h_eegchan = get(h_eeg, 'children');
set(h_eegchan(eeg_chan_num(selected_chan)),'color', 'r')

axes(h_topo);
hold on
z = 2.1; %emarkerheight

plot3(x(eeg_chan_num(selected_chan)), y(eeg_chan_num(selected_chan)),...
    ones(size(x(eeg_chan_num(selected_chan))))*z,...
    '.','color','r','markersize', 20, 'linewidth', 1);

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn, fig, h_eeg, chan_sel_status, eeg_chan_num, sel_chan_names, liststring, h_topo, x,y) %#ok

% Highlight channel in eegplot
h_eegchan = get(h_eeg, 'children');

listval = get(listbox,'value');

for i = 1:length(listval)
    selected_chan(i) = find(strcmp(liststring{listval(i)}, sel_chan_names) == 1);
end

if chan_sel_status(selected_chan) == 0
    chan_sel_status(selected_chan) = 1;
else
    chan_sel_status(selected_chan) = 0;
end

axes(h_eeg); % shift focus to eegplot

set(h_eegchan(eeg_chan_num(find(chan_sel_status == 0))),'color', [0 0 0.4])
set(h_eegchan(eeg_chan_num(find(chan_sel_status == 1))),'color', 'r')

% set(h_eegchan(find(chan_sel_status == 0)),'color', 'b')
% set(h_eegchan(find(chan_sel_status == 1)),'color', 'r')

% Highlight channel loc in topoplot
axes(h_topo); % shift focus to topoplot
hold on
z = 2.1; %emarkerheight

sel_locs = findobj(gca,'color','r');
set(sel_locs,'color','k','markersize',10);

plot3(x(eeg_chan_num(find(chan_sel_status == 1))), y(eeg_chan_num(find(chan_sel_status == 1))),...
    ones(size(x(eeg_chan_num(find(chan_sel_status == 1)))))*z,...
    '.','color','r','markersize', 20, 'linewidth', 1);

% plot3(x(find(chan_sel_status == 1)), y(find(chan_sel_status == 1)),...
%     ones(size(x(find(chan_sel_status == 1))))*z,...
%     '.','color','r','markersize', 20, 'linewidth', 1);

% return focus to listbox
figure(fig)

% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);       
else
    if length(get(listbox,'string'))==length(get(listbox,'value'))
        set(selectall_btn,'enable','off')
    else
        set(selectall_btn,'enable','on')
    end    
end
