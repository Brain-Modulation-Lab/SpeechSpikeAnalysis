function [com,comb,f,R,zg]=sliding_fft(signal, fs,Events, baseline,stat,frange,tstart,tend)
%function [com,comb,f,R,zg]=sliding_fft(signal, fs,Events, baseline,stat,frange,tstart,tend)

if size(signal,1)>size(signal,2)
    signal=signal';
end
ch=size(signal,1);

winlen=256; winshift=floor(winlen*0.9); nfft=winlen*2;
f= linspace(0,fs/2,floor(nfft/2)+1);
fi_lb = find(f >= frange(1),1,'first'); if isempty(fi_lb) fi_lb=1; end;
fi_ub = find(f <= frange(2),1,'last'); if isempty(fi_ub) fi_ub = length(f); end;
f=f(fi_lb:fi_ub);
hwin=hann(winlen);
prestim=winlen*10; poststim=winlen*10; Events=round(Events*fs);
baseline=round(baseline*fs); baseline_dur=winlen*5;
wstart=1:winshift:prestim+poststim+1-winlen;
signal=squeeze(num2cell(signal',1));
wind=zeros(winlen,length(wstart));
wind(:,1)=1:winlen;
for i=2:length(wstart)
    wind(:,i)=wstart(i)+1:wstart(i)+winlen;
end
bind=[];
for b=1:length(baseline)
    bind=cat(2,bind,baseline(b)-baseline_dur:baseline(b));
end
tind=[];
for t=1:length(Events)
    tind=cat(2,tind,Events(t)-prestim:Events(t)+poststim);
end

dx=length(tind)/length(Events);   dz=length(Events);
tr=cellfun(@(x) reshape(x(tind,:)',dx,dz),signal,'UniformOutput', false);
dyb=length(bind)/length(baseline);   dzb=length(baseline);
base=cellfun(@(x) reshape(x(bind,:)',dyb,dzb),signal,'UniformOutput', false);

com=zeros(length(f),length(wstart),dz,ch);
for i=1:length(tr)
    x=tr{i};
    tmp=zeros(nfft,length(wstart),dz);
    for i2=1:size(x,2)
        y=x(:,i2);
        tmp(:,:,i2)=fft(bsxfun(@times,y(wind),hwin),nfft)/winlen;
    end
    com(:,:,:,i)=tmp(frange(1):frange(2),:,:);
end
wstartb=1:winshift:baseline_dur+1-winlen;
windb=zeros(winlen,length(wstartb));
windb(:,1)=1:winlen;
for i=2:length(wstartb)
    windb(:,i)=wstartb(i)+1:wstartb(i)+winlen;
end
comb=zeros(length(f),length(wstartb),dzb,ch);
for i=1:ch
    x=base{i}(1:end,:);
    tmp=zeros(nfft,length(wstartb),dzb);
    for i2=1:size(x,2)
        y=x(:,i2);
        tmp(:,:,i2)=fft(bsxfun(@times,y(windb),hwin),nfft)/winlen;
    end
    comb(:,:,:,i)=tmp(frange(1):frange(2),:,:);
end



parfor_progress(stat.surrn);
parfor s=1:stat.surrn
%surrogates
sind=[];
sevents=randperm(round(tend-tstart-2*poststim),dz)+round(tstart);
sevents(sevents<=prestim)=sevents(sevents<=prestim)+prestim;
for t=1:length(sevents)
    sind=cat(2,sind,sevents(t)-prestim:sevents(t)+poststim);
end
% sbaseline=sevents-2.5*fs;
% bind=[];
% for b=1:length(sbaseline)
%     bind=cat(2,bind,sbaseline(b)-baseline_dur:sbaseline(b));
% end
% sbase=cellfun(@(x) reshape(x(bind,:)',dxb,dz),signal,'UniformOutput', false);

sr=cellfun(@(x) reshape(x(sind,:)',dx,dz),signal,'UniformOutput', false);
[coms(:,:,:,s)]=whysorandom(sr,hwin,wstart,dz,nfft,winlen,wind,frange);
parfor_progress;
end
parfor_progress(0);

parfor i=1:ch
[R{i}, zg{i}]=cluster_sig(bsxfun(@minus,mean(abs(com(:,:,:,i)),3),mean(mean(abs(comb(:,:,:,i)),2),3)),bsxfun(@minus,squeeze(coms(:,:,i,:)),mean(mean(abs(comb(:,:,:,i)),2),3)),stat,2);
end
end

function [coms]=whysorandom(sr,hwin,wstart,dz,nfft,winlen,wind,frange)



for i=1:length(sr)
    x=sr{i};
    tmp=zeros(nfft,length(wstart),dz);
    for i2=1:size(x,2)
        y=x(:,i2);
        tmp(:,:,i2)=fft(bsxfun(@times,y(wind),hwin),nfft)/winlen;
    end
    coms(:,:,i)=mean(abs(tmp(frange(1):frange(2),:,:)),3);
end


end