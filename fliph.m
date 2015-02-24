function N=fliph(M)
[nv nh]=size(M);
for k=1:nv
for l=1:nh
N(k,l)=M(nv-k+1,l);
end
end