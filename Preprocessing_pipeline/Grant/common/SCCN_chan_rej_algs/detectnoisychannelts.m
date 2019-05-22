function [badchanindts,badchanlogts,c_cellts,IND,maxcts]=detectnoisychannelts(EEG)
% 
warning('OFF', 'MATLAB:divideByZero');
szdata=size(EEG.data(:,:));
pos=[cell2mat({EEG.chanlocs(1:szdata(1)).X}); cell2mat({EEG.chanlocs(1:szdata(1)).Y}); cell2mat({EEG.chanlocs(1:szdata(1)).Z});];
ind=[1:szdata(1)];
for n=1:szdata(1)
    npos=pos(:,n);
    Rnpos=repmat(npos,[1 size(pos,2)]);
    dist=sqrt(sum((Rnpos-pos).^2));
    i=find(dist<.25 & dist>0); % was .5 before here the distance between electrodes is defined to look for correlated data; .1 = 10 cm
    IND{n}=i;
end
numt=floor(szdata(2)/256);
badchanlogts=zeros(szdata(1),numt);

EEG.data = EEG.data(:,:); % nima

for t=1:numt
    %c_vec=[];
  %  t
    data=EEG.data(:,[(256*t -255):256*t]);
%    fprintf('Calculating correlation matrix ...\n');
    C=corrcoef([data(:,:)']);
 %   fprintf('done.\n');
    for n=1:szdata(1)
        i=IND{n};
        
%         if isempty(i)
%            fprintf('nearby channels are empty at %d\n',n)
%         end;
        
        %C=corrcoef([data(n,:)' data(i,:)']);
        %c=C(1,2:end);
        c= C(n,i);
        c_cell{n}=c;
        %c_vec=[c_vec c];       
       % n
       if ~isempty(c)
           maxc(n)=max(c);
       end;
    end
    badchanind=find(abs(maxc)<.4);
    badchanindts{t}=badchanind;
    maxcts{t}=maxc;
    c_cellts{t}=c_cell;
    badchanlogts(badchanind,t)=1;
end
warning('ON', 'MATLAB:divideByZero');