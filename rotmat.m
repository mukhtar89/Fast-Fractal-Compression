function N=rotmat(M)
[nv nh]=size(M);
for k=1:nv
for l=1:nh
N(k,l)=M(l,nv-k+1);
end
end