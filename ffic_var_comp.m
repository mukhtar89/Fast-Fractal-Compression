% Fast Fractal compression using Variance Ordering Method

clc;
clf;
% Set timers
begrun=clock
cpu=cputime
M=imread('lena_gray_256.tif');
[sv sh]=size(M);
if sv~=sh
    display('Matrix is not square');
    return
end

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

% Find variances of domain block for different scalings
Vd=zeros(nd,nd,4);
for i=1:nd
    i1=(i-1)*rsize+1;
    i2=i*rsize;
    for j=1:nd
        j1=(j-1)*rsize+1;
        j2=j*rsize;
        D=M1(i1:i2,j1:j2);
        D=double(D);
        Z=mean(mean(D));
        Z=D-Z;
        for n=1:4
            Vd(i,j,n)=(s(n)^2)*sum(sum(Z.^2));
        end
    end
end



% Sorting variances of domain block
Vd_ord=zeros(4,nd*nd*4);
Vd_ord=Vd_ord';
w=Vd(:);
Vd_ord(:,1)=w;
e=1:nd;
e=e';
e=e*ones(1,nd);
e=e(:);
e=e*ones(1,4);
e=e(:);
f=1:nd;
f=f'*ones(1,nd);
f=f';
f=f(:);
f=f*ones(1,4);
f=f(:);
g=[1 2 3 4];
g=g';
g=g*ones(1,nd*nd);
g=g';
g=g(:);
Vd_ord(:,2)=e;
Vd_ord(:,3)=f;
Vd_ord(:,4)=g;


for k=1:(nd*nd)
    for i=1:(nd*nd-1)
        temp=zeros(1,4);
        if Vd_ord(i,1)<Vd_ord(i+1,1)
            for j=1:4
                temp(j)=Vd_ord(i,j);
                Vd_ord(i,j)=Vd_ord(i+1,j);
                Vd_ord(i+1,j)=temp(j);
            end
        end
    end
end


% Begin batch runs
clear T;

% Compare the range blocks and scaled domain blocks.
% k,l - used to cycle through blocks Rkl.
for k=1:nr
    k1=(k-1)*rsize+1;
    k2=k*rsize;
    for l=1:nr
        l1=(l-1)*rsize+1;
        l2=l*rsize;
        R=M(k1:k2,l1:l2);
        R=double(R);
        % Offset o is the average in the block Rkl
        O=mean(mean(R));
        Rm=R-O;
        Vr=sum(sum(Rm.^2));
        min=10000;
        u=1;
        for i4=1:(nd*nd)
            er=(sqrt(Vr)-sqrt(Vd_ord(i4,1)))^2;
            if er<min
                min=sqrt(er);
                u=i4;
            end
        end
        
        % Initialize error to large value
        i0=0;
        j0=0;
        m0=0;
        g0=0;
        s0=0;
        
        dmin=10^9;
        stf=0;
        q=zeros(1,2);
        q(1)=u-1;
        q(2)=u;
        while stf==0
            dir=1;
            while (dir>=0)&(stf==0)
                if ((sqrt(Vd_ord(q(dir+1),1))-sqrt(Vr))^2)>=dmin^2
                    stf=2;
                else
                    i=Vd_ord(q(dir+1),2);
                    j=Vd_ord(q(dir+1),3);
                    n=Vd_ord(q(dir+1),4);
                    
                    % Now cycle through each Domain Dij
                    i1=(i-1)*rsize+1;
                    i2=i*rsize;
                    j1=(j-1)*rsize+1;
                    j2=j*rsize;
                    
                    
                    bigM=zeros(rsize,rsize,8);
                    D=M1(i1:i2,j1:j2);
                    del_g=O-s(n)*mean(mean(D));
                    D=D+del_g;
                    bigM(:,:,1)=D;
                    tmp=rotmat(D);
                    bigM(:,:,2)=tmp;
                    tmp=rotmat(tmp);
                    bigM(:,:,3)=tmp;
                    tmp=rotmat(tmp);
                    bigM(:,:,4)=tmp;
                    bigM(:,:,5)=fliph(D);
                    bigM(:,:,6)=flipv(D);
                    bigM(:,:,7)=D';
                    bigM(:,:,8)=rotmat(rotmat(D'));
                    
                    % Test each transformation
                    for m=1:8
                        D=bigM(:,:,m);
                        sum_dist=sum(sum((R-D).^2));
                        dist=sqrt(sum_dist);
                        if dist<dmin
                            dmin=dist;
                            i0=i;
                            j0=j;
                            m0=m;
                            s0=s(n);
                            g0=del_g;
                        end
                    end
                end
                q(dir+1)=q(dir+1)+((-1)^dir);
                if ((q(dir+1)<1)|(q(dir+1)>(nd*nd*4)))==1
                    stf=2;
                end
            end
        end
        T(k,l,:)=[i0 j0 m0 s0 g0];
    end
end

% Stop the clock, store computation time in tim
% and elapsed cpu time in cpu0.
cpu0=cputime-cpu
stoprun=clock
tim=etime(begrun,stoprun)
% Save data in mat file - need to change the name after each use.
save 'gs_ffic' sv rsize T tim cpu0;