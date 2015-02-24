clc;
clf;
M=imread('grayscale.gif');
[sv sh]=size(M)
if sv~=sh
    display('Matrix is not square');
    return
end

t1=[0 0 1 7 5 5 0 1 0 0;
    0 1 0 5 7 0 1 0 1 0;
    0 5 7 0 2 2 0 0 5 0;
    0 7 5 2 0 0 7 5 0 0;
    0 1 0 5 0 0 4 5 7 0;
    0 0 1 0 2 4 0 7 5 0;
    0 1 0 0 2 7 5 0 3 0;
    0 0 1 7 0 5 7 3 0 0;
    0 0 0 0 0 0 0 0 0 0];

% Begin batch runs
for irn=1:10
    clear T;
    % Set timers
    begrun=clock;
    cpu=cputime;
    min0=10*irn;
    rsize=4;
    nd=sv/(rsize*2);
    nr=sv/rsize;
    % Scale the Domain Blocks
        for i=0:rsize*nd-1
            for j=0:rsize*nd-1
                M1(i+1,j+1)=(M(2*i+1,2*j+1)+M(2*i+2,2*j+1)+M(2*i+1,2*j+1)+M(2*i+2,2*j+2))/4;
            end
        end

    % Create monster matrix containing all possible 2D transformations
    % of the domain blocks. Store in multidimensional matrix bigM.
    for i=1:nd
        i1=(i-1)*rsize+1;
        i2=i*rsize;
        for j=1:nd
            j1=(j-1)*rsize+1;
            j2=j*rsize;
            D=M1(i1:i2,j1:j2);
            %D=D-mean(mean(D));
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
    for k=1:nr
        k1=(k-1)*rsize+1;
        k2=k*rsize;
        for l=1:nr
            l1=(l-1)*rsize+1;
            l2=l*rsize;
            R=M(k1:k2,l1:l2);
            % Offset o is the average in the block Rkl
            R=double(R);
            Rtype=-1;
            thres=5;
            
            diagrdn=zeros(4,4);
            diagrup=zeros(4,4);
            for e=1:4
                for f=1:4
                    if(f<=e)
                        diagrdn(e,f)=R(e,f);
                    else
                        diagrup(e,f)=R(e,f);
                    end
                end
            end
            diagldn=zeros(4,4);
            diaglup=zeros(4,4);
            for e=1:4
                for f=1:4
                    if(f<=(4-e))
                        diaglup(e,f)=R(e,f);
                    else
                        diagldn(e,f)=R(e,f);
                    end
                end
            end
            
            
            if((abs(sum(sum(R(:,1:2)))-sum(sum(R(:,3:4))))<thres)&&(abs(sum(sum(R(1:2,:)))-sum(sum(R(3:4,:)))))<thres)
                Rtype=0
            elseif sum(sum(R(:,1:2)))<sum(sum(R(:,3:4)))
                Rtype=1
            elseif sum(sum(R(:,1:2)))>sum(sum(R(:,3:4)))
                Rtype=2
            elseif sum(sum(R(1:2,:)))<sum(sum(R(3:4,:)))
                Rtype=3
            elseif sum(sum(R(1:2,:)))>sum(sum(R(3:4,:)))
                Rtype=4
            elseif sum(sum(diaglup))>sum(sum(diagldn))
                Rtype=5
            elseif sum(sum(diaglup))<sum(sum(diagldn))
                Rtype=6
            elseif sum(sum(diagrup))<sum(sum(diagrdn))
                Rtype=7
            elseif sum(sum(diagrup))>sum(sum(diagrdn))
                Rtype=8
            else
                Rtype=9
            end
        
            % Initialize error to large value
            minerr=1000;
            i0=0;
            j0=0;
            m0=0;
            t0=0;
            
                % Now cycle through each Domain Dij
            if minerr>min0
            for i=1:nd
                    if minerr>min0
                        i1=(i-1)*rsize+1;
                        i2=i*rsize;
                        for j=1:nd
                            if minerr>min0
                                j1=(j-1)*rsize+1;
                                j2=j*rsize;
                                
                                Dtype=-1;
                                % Test each transformation
                                for m=1:8
                                    if minerr>min0
                                        D=bigM(i1:i2,j1:j2,m);
                                        
                                        D1=dct(D);
                                        E=sum(sum(D1^2))/(rsize^2);
                                        D2=D1/E;
                                        
                                        if(D2(1,1)>0.9)
                                            Dtype=0
                                        elseif (abs(D1(1,2))-abs(D1(2,1))>32)&(D1(1,2)<0)
                                                Dtype=1
                                        elseif (abs(D1(1,2))-abs(D1(2,1))>32)&(D1(1,2)>0)
                                                Dtype=2
                                        elseif (abs(D1(2,1))-abs(D1(1,2))>32)&(D1(1,2)<0)
                                                Dtype=3
                                        elseif (abs(D1(2,1))-abs(D1(1,2))>32)&(D1(1,2)>0)
                                                Dtype=4
                                        elseif (abs(abs(D1(1,2))-abs(D1(2,1)))<24)&(D1(1,2)>0)&(D1(2,1)>0)
                                                Dtype=5
                                        elseif (abs(abs(D1(1,2))-abs(D1(2,1)))<24)&(D1(1,2)<0)&(D1(2,1)<0)
                                                Dtype=6
                                        elseif (abs(abs(D1(1,2))-abs(D1(2,1)))<24)&(D1(1,2)>0)&(D1(2,1)<0)
                                                Dtype=7
                                        elseif (abs(abs(D1(1,2))-abs(D1(2,1)))<24)&(D1(1,2)<0)&(D1(2,1)>0)
                                                Dtype=8
                                        else
                                            Dtype=9
                                        end
                    
                                        R_avg=mean(mean(R));
                                        D_avg=mean(mean(D));
                                        if(Rtype~=0)
                                            t=t1(Rtype,Dtype+1)
                                        else
                                            t=R_avg
                                        end
                                        
                                        a_num=0;
                                        a_den=0;
                                        for e=1:4
                                            for f=1:4
                                                a_num=a_num+((R(e,f)-R_avg)*(D(e,f)-D_avg));
                                                a_den=a_den+((D(e,f)-D_avg)^2);
                                            end
                                        end
                                        a=a_num/a_den;
                                        b=R_avg-a*D_avg;
                                        
                                        er=0;
                                        for e=1:4
                                            for f=1:4
                                                er=er+R(e,f)-a*D(e,f)-b;
                                            end
                                        end
                                        if er<minerr
                                                    minerr=er
                                                    i0=i;
                                                    j0=j;
                                                    m0=m;
                                                    a0=a;
                                                    b0=b;
                                                    t0=t;
                                        end
                                    end
                                end
                                                    end
                                                end
                                             end
                        end
            end
                         T(k,l,:)=[i0 j0 m0 a0 b0 t0];
                      end   
            end
            
    % Stop the clock, store computation time in tim
    % and elapsed cpu time in cpu0.
    cpu0=cputime-cpu;
    stoprun=clock;
    tim=etime(begrun,stoprun);
    % Save data in mat file - need to change the name after each use.
    switch irn
    case 1,
    save 'gs4_1' sv rsize T tim cpu0;
    case 2,
    save 'gs4_2' sv rsize T tim cpu0;
    case 3,
    save 'gs4_3' sv rsize T tim cpu0;
    case 4,
    save 'gs4_4' sv rsize T tim cpu0;
    case 5,
    save 'gs4_5' sv rsize T tim cpu0;
    case 6,
    save 'gs4_6' sv rsize T tim cpu0;
    case 7,
    save 'gs4_7' sv rsize T tim cpu0;
    case 8,
    save 'gs4_8' sv rsize T tim cpu0;
    case 9,
    save 'gs4_9' sv rsize T tim cpu0;
    case 10,
    save 'gs4_10' sv rsize T tim cpu0;
    end
end