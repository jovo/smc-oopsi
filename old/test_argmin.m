a=rand(10,1);
len=1e5;
b=zeros(len,1);
c=zeros(len,1);
d=zeros(len,1);
for i=1:len
    b(i)=argmin(a);
    c(i)=find(a,1);
    [foo, d(i)]=min(a(:));
end