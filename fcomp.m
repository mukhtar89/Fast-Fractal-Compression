clc;
clf;
% Set timers
begrun=clock
cpu=cputime
M=imread('grayscale.gif');
[sv sh]=size(M);
if sv~=sh
    display('Matrix is not square');
    return
end
% Begin batch runs


min0=100;
rsize=4;
nd=sv/rsize/2;
nr=sv/rsize;
% Scale the Domain Blocks
for i=1:rsize*nd
    for j=1:rsize*nd
        M1(i,j)=mean(mean(M((i-1)*2+1:i*2,(j-1)*2+1:j*2)));
    end
end
% Matrix of 4 possible scalings to transform grayscale
s=[0.45 0.60 0.80 1];
% Create monster matrix containing all possible 2D transformations
% of the domain blocks. Store in multidimensional matrix bigM.
for i=1:nd
    i1=(i-1)*rsize+1;
    i2=i*rsize;
    for j=1:nd
        j1=(j-1)*rsize+1;
        j2=j*rsize;
        D=M1(i1:i2,j1:j2);
        D=D-mean(mean(D));
        bigM(i1:i2,j1:j2,1)=D;
        tmp=rotmat(D);
        bigM(i1:i2,j1:j2,2)=tmp;
        tmp=rotmat(tmp);
        bigM(i1:i2,j1:j2,3)=tmp;
        tmp=rotmat(tmp);
        bigM(i1:i2,j1:j2,4)=tmp;
        bigM(i1:i2,j1:j2,5)=fliph(D);
        bigM(i1:i2,j1:j2,6)=flipv(D);
        bigM(i1:i2,j1:j2,7)=D';
        bigM(i1:i2,j1:j2,8)=rotmat(rotmat(D'));
    end
end
% Compare the range blocks and scaled domain blocks.
% k,l - used to cycle through blocks Rkl.
clear T;
for k=1:nr
    k1=(k-1)*rsize+1;
    k2=k*rsize;
    for l=1:nr
        l1=(l-1)*rsize+1;
        l2=l*rsize;
        R=M(k1:k2,l1:l2);
        % Offset o is the average in the block Rkl
        o=mean(mean(R));
        R=R-o;
        R=double(R);
        % Initialize error to large value
        minerr=10000;
        i0=0;
        j0=0;
        m0=0;
        if minerr>min0
            % Now cycle through each Domain Dij
            for i=1:nd
                if minerr>min0
                    i1=(i-1)*rsize+1;
                    i2=i*rsize;
                    for j=1:nd
                        if minerr>min0
                            j1=(j-1)*rsize+1;
                            j2=j*rsize;
                            % Test each transformation
                            for m=1:8
                                if minerr>min0
                                    D=bigM(i1:i2,j1:j2,m);
                                    D=double(D);
                                    % Try the four gray scalings
                                    for n=1:4
                                        if norm(s(n)*D-R)<minerr
                                            minerr=norm(s(n)*D-R);
                                            i0=i;
                                            j0=j;
                                            m0=m;
                                            z=s(n);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        T(k,l,:)=[i0 j0 m0 z o];
    end
end
% Stop the clock, store computation time in tim
% and elapsed cpu time in cpu0.
cpu0=cputime-cpu
stoprun=clock
tim=etime(begrun,stoprun)
% Save data in mat file - need to change the name after each use.
save 'gs_norm' sv rsize T tim cpu0;
