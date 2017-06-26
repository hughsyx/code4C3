function [CAC1,CAC2,CAC] = C3_SC(src,rcv,lat,lon,file_name,v_src,v_rcv)


tic
npt1 = 3000;   % length of the coda wave
vel = 2.8;

% get the index of the source and receiver in file_name
for i=1:length(file_name)
    if strcmp(file_name(i).name,src)
        ind_src = i;
    end
    if strcmp(file_name(i).name,rcv)
        ind_rcv = i;
    end
end


Npt = round(size(v_rcv,1)+1)/2;
npt = 8000; % length of the C1 that going to use (on each side (causal and anti-causal))
dur = [-npt+1:npt-1]+Npt;
[b,a] = butter(4,[(2*1/10*0.2),(2*1/3*0.2)],'bandpass');
V_src = filtfilt(b,a,v_src(dur,:));
V_rcv = filtfilt(b,a,v_rcv(dur,:));


CAC1=0;CAC2=0;CAC=0;
dist = distance([lat(ind_src),lon(ind_src)],[lat(ind_rcv),lon(ind_rcv)])*111.11; % distance between the virtual source and receiver


for i=1:length(file_name)
    CBA = V_src(:,i); CBC = V_rcv(:,i);
    
        dist1 = distance([lat(ind_src),lon(ind_src)],[lat(i),lon(i)])*111.11;
        dist2 = distance([lat(ind_rcv),lon(ind_rcv)],[lat(i),lon(i)])*111.11;
        if max(dist1,dist2) > 450
            continue
        end

        CBA1=fliplr(CBA(1:npt)')';  % separate anti-causal part and causal part
        CBC1=fliplr(CBC(1:npt)')';
        CBA2=CBA(npt:end);
        CBC2=CBC(npt:end);
        
        [CBA_coda1,CBA_coda2,CBC_coda1,CBC_coda2] = get_coda(CBA1,CBA2,CBC1,CBC2,dist1,dist2,npt1,vel); % get coda part of C1

        if max(abs(CBA_coda1))==0 | max(abs(CBC_coda1))==0
            continue
        end
        
        %get c3 of the cacusal and anti-caculsal separately
        cc1 = CANIR2(CBC_coda1,CBA_coda1,10);
        CAC1 = CAC1+cc1;
        cc2 = CANIR2(CBC_coda2,CBA_coda2,10);
        CAC2 = CAC2+cc2;
end
toc