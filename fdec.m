
load 'gs_norm'
% Initialize matrix
M=100*ones(sv);
% Start Iteration
for iter=1:10
    % Enter range block size used in fcomp
    rsize=4;
    nd=sv/rsize/2;
    nr=sv/rsize;
    % Rescale Domain Blocks
    for i=1:rsize*nd
        for j=1:rsize*nd
            M1(i,j)=mean(mean(M((i-1)*2+1:i*2,(j-1)*2+1:j*2)));
        end
    end
    % Transform Domain Block Using T matrix
    for k=1:nr
        k1=(k-1)*rsize+1;
        k2=k*rsize;
        for l=1:nr
            l1=(l-1)*rsize+1;
            l2=l*rsize;
            i0 = T(k,l,1);
            j0 = T(k,l,2);
            m0 = T(k,l,3);
            s0 = T(k,l,4);
            g0 = T(k,l,5);
            i1 = (i0-1)*rsize+1;
            i2 = i0*rsize;
            j1 = (j0-1)*rsize+1;
            j2 = j0*rsize;
            D = M1(i1:i2,j1:j2);
            if m0==2
                D=rotmat(D);
            elseif m0==3
                D=rotmat(rotmat(D));
            elseif m0==4
                D=rotmat(rotmat(rotmat(D)));
            elseif m0==5
                D=fliph(D);
            elseif m0==6
                D=flipv(D);
            elseif m0==7
                D=D';
            elseif m0==8
                D=rotmat(rotmat(D'));
            end
            R=s0*D+g0*ones(size(D));
            MM(k1:k2,l1:l2)=R;
        end
    end
    M=MM;
end
% Output Image which is in M
imagesc(M)
colormap(gray);