function ANIRW=CANIR2(vb,va,npt)
% calculate the ANIR of each window pair
% va is the virtual source
% npt is the number used for running window average

len1=length(va);len2=length(vb);
len=len1+len2-1;
fva=fft(va,len);
fvb=fft(vb,len);

eps=1e-9;

fva2=RWA(fva,npt); fvb2=RWA(fvb,npt);
p1=ones(size(fva))*eps; p2=ones(size(fvb))*eps;
for i=1:len
    p1(i)=sum(fva2(i:i+npt-1))/npt;
    p2(i)=sum(fvb2(i:i+npt-1))/npt;
end

ANIRWF=fvb.*conj(fva)./(p1.*p1);
ANIRW=real(fftshift(ifft(ANIRWF)));

function v2 = RWA(v,npt)
    temp = zeros(1,npt-1);
    len = length(v);
    v2 = zeros(len+npt-1,1);
    v2(1:len) = v;
    temp=v(end-npt+1:end-1);
    v2(len+1:end)=fliplr(temp);
    v2=abs(v2);
end
end
