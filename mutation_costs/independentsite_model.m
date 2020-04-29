%% Independent-site model

h=zeros(402,1);
DE_av=zeros(402,1);
DE_calculated=zeros(402,1);

Mend=130000;
sb = 360/408;

M = size(X,1);
f = (1/M)*(sum(X));

%% Field of the standard indepedent-site model
hi = real(atanh(2*f-1));

for i = 1:402
    DE_calculated(i) = (2*f(i)-1)*hi(i);       % estimation conservation
end

%% Field with background frequency
for i = 1:402
    h(i)=log(f(i)*(1-sb)/(sb*(1-f(i)))); %field
end

for i=1:402
 DE_av(i) = (2*f(i)*(1-sb)/((1-f(i))*sb + f(i)*(1-sb)) - 1)*h(i);
end

%% Conservation of Halabi et al.
cons_Halabi=zeros(402,1);
for i=1:402
 cons_Halabi(i) = f(i)*log(f(i)/sb) + (1-f(i))*log((1-f(i))/(1-sb));
end




