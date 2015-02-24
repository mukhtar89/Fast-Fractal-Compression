function N=flipv(M)
[nv nh]=size(M);
for k=1:nv
for l=1:nh
N(k,l)=M(k,nv-l+1);
end
end