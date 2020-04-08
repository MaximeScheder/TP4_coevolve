%% Independent-site model

h=zeros(402,1);
DE_tmp=zeros(402,1);
DE_calculated=zeros(402,1);
Mend=130000;
sb = 360/408;

M = size(X,1);
f = (1/M)*(sum(X));

%hi = real(atanh(2*f-1));

for m = 1:Mend
for i = 1:402
    h(i)=log(f(i)*(1-sb)/(sb*(1-f(i)))); %field
    DE_tmp(i) = (2*X(m,i)-1)*h(i);       %mutation cost
end
DE_calculated = DE_calculated + DE_tmp;
end
DE_calculated = DE_calculated/Mend;
