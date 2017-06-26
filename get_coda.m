function  [Coda1,Coda2,Coda3,Coda4]=get_coda(C1,C2,C3,C4,dist1,dist2,npt,vel)
% get coda part of C1
%
    c=vel;
    f = 5; % sampling rate of the input
    ind1=round(dist1/c*2*f); ind2=round(dist2/c*2*f);
    % truncate the signal start from 2 times travel time
    ind = max(ind1,ind2);

    Coda1=C1(ind:ind+npt-1);
    Coda2=C2(ind:ind+npt-1);
    Coda3=C3(ind:ind+npt-1);
    Coda4=C4(ind:ind+npt-1);
end
