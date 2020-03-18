function FitnessCOOPPB(Minit, Mend, filename, saveName)

%PERIODIC BOUNDARIES

 config=load(filename);
% config=load('/Users/maximescheder/Documents/EPFL/MA1_1/Labo/DCA_ex/shear05_402_PB.dat');
% config=load('DCA_inf.mat');
% config=config.samples_DCA;

%config=ones(1000,408);
lim=Mend; %was at 5
configfull=[config(1:lim,1:12), ones(lim,1), config(1:lim,13:14), ones(lim,1), config(1:lim,15:16), ones(lim,1), config(1:lim,17:397), ones(lim,3), config(1:lim,398:402)];

consensus=load('C:\Users\msche\OneDrive\Documents\EPFL\TP4\DCA_ex\consensus_shearPB.dat');
    
xl=12;
yl=12;
nsp=360;
beta_ini=0.0001;

allos = 0;  % 0 for cooperative, 1 for geometric
dd    = 1;
dt    = 1.;
kweak = 0.0001;
drbin = 0.2;
topnb = 36;

lmem  = 32;
%dmem  = 'uint32';

np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
temprange = 1;
eps   = 0.2;   % distort size
%confint = 100;
%slow = 10;
warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl;
dim   = 2;
xbound = xl;
ybound = yl*sqrt(3)/2;
nmax  = 3*num-2*xl;
nflip = ones(1,temprange);
nf    = ones(1,temprange);
beta  = 1/beta_ini;
for tt = 1:temprange
    nflip(tt) = temprange+1-tt;
end

%% define perturb and target displacements
id  = (xl-np)/2+(1:np);     % indices of the displaced particles
idd = dim*id(1)-1:dim*id(end);  %
it  = (xl-nt)/2+(yl-1)*xl+(1:nt);  % indices of the target particles
itt = dim*it(1)-1:dim*it(end);
idt = [idd,itt];
nid = setxor(1:dim*num,idd);
nit = setxor(1:dim*num,itt);
nidt= setxor(1:dim*num,[idd,itt]);
ndd = 2*dim+1:dim*num;
pts  = zeros(dim*num,dim);  % translation on perturbing nodes
tts  = zeros(dim*num,dim);  % translation on target nodes
pert = zeros(dim*num,1);  % perturb displacement
targ = zeros(dim*num,1);  % target displacement
ptx  = zeros(dim*num,1);  % translation in x direction on perturbing nodes
pty  = zeros(dim*num,1);  % translation in y direction on perturbing
prot = zeros(dim*num,1);  % rotation on perturbing
ttx  = zeros(dim*num,1);  % translation in x direction on target nodes
tty  = zeros(dim*num,1);  % translation in y direction on target
trot = zeros(dim*num,1);  % rotation on target
tx   = zeros(dim*num,1);  % translation in x direction on target nodes
ty   = zeros(dim*num,1);  % translation in y direction on target
rot  = zeros(dim*num,1);  % rotation on target
fixs = zeros(1,(np+nt)*2);  % labels of springs to be fixed
cf   = 0;

%% construct the embedded network
bondnb = zeros(3*num,2);
nnb    = 0;
bondwk = zeros(6*num,2);  % to all next neareat neighbors
nwk    = 0;
posx = zeros(1,num);
posy = zeros(1,num);
for nn=1:num
    ny = ceil(nn/xl);
    nx = nn - (ny-1)*xl;
    n = nn;
    %nx = mod(nx+7,xl)+1;
    %n  = (ny-1)*xl+nx;
    posx(n) = nx-1+mod(ny-1,2)/2;
    posy(n) = sqrt(3)/2*ny - sqrt(3)/4;
    %
    if mod(ny,2)
        if mod(nx+mod(floor(ny/2),2),2)
            posx(n) = posx(n);
            posy(n) = posy(n);
        else
            posx(n) = posx(n) + eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    else
        if mod(nx+mod(ny/2-1,2),2)
            posx(n) = posx(n);
            posy(n) = posy(n) - eps;
        else
            posx(n) = posx(n) - eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    end
    %}
    %theta = rand(1);
    %posx(n) = posx(n) + eps*cos(theta);
    %posy(n) = posy(n) + eps*sin(theta);
    %if nx<xl
    nnb = nnb+1;
    bondnb(nnb,1) = n;
    bondnb(nnb,2) = n+1-xl*(nx==xl);  % periodic in horizontal direction
    %end
    % fix the connections between the perturbing nodes and target nodes
    if ismember(n,id)&&ismember(n+1-xl*(nx==xl),id)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ismember(n,it)&&ismember(n+1-xl*(nx==xl),it)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ny<yl  % the top line is not periodic
        %if nx<xl ||mod(ny,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl+(1-xl*(nx==xl))*mod(ny-1,2); %-num*(ny==linnum)
        %end
        %if nx>1 ||mod(ny+1,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl-(1-xl*(nx==1))*mod(ny,2);  %-num*(ny==linnum)
        %end
    end
    % weak bonds connect to the next nearest and next next nearest
    % next nearest neighbors
    if ny<yl
        %if nx<xl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl+1-xl*(nx==xl)+(1-xl*(nx==xl-1))*mod(ny-1,2);  %-num*(ny==linnum)
        %end
        %if nx>2
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl-1+xl*(nx==1)-(1-xl*(nx==2))*mod(ny,2);  %-num*(ny==linnum)
        %end
    end
    if ny<yl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl;
        % next-next nearest neighbors
        %if nx<xl
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl+1-xl*(nx==xl);
        %end
        %if nx>1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl-1+xl*(nx==1);
        %end
    end
   %if nx<xl-1
    nwk = nwk+1;
    bondwk(nwk,1) = n;
    bondwk(nwk,2) = n+2-xl*(nx>=xl-1);
   % end
end
bondnb = bondnb(1:nnb,:);
bondwk = bondwk(1:nwk,:);
fixs   = fixs(1:cf);

%% initialize the perturbing and target displacements
pert(2*id(1)) = 1*dd;
targ(2*it(1)) = -1*dt;
for i = 2:np
    dx  = posx(id(i))-posx(id(i-1));
    dy  = posy(id(i))-posy(id(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    pert(2*id(i)-1) = pert(2*id(i-1)-1)*(dx^2-dy^2)+2*pert(2*id(i-1))*dx*dy;
    pert(2*id(i))   = pert(2*id(i-1))*(dy^2-dx^2)+2*pert(2*id(i-1)-1)*dx*dy;
end
for i = 2:nt
    dx  = posx(it(i))-posx(it(i-1));
    dy  = posy(it(i))-posy(it(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    targ(2*it(i)-1) = targ(2*it(i-1)-1)*(dx^2-dy^2)+2*targ(2*it(i-1))*dx*dy;
    targ(2*it(i))   = targ(2*it(i-1))*(dy^2-dx^2)+2*targ(2*it(i-1)-1)*dx*dy;
end
% translational and rotational parts at the stimulus site and target site
pcen = [mean(posx(id)),mean(posy(id))]; % center of perturbing nodes
tcen = [mean(posx(it)),mean(posy(it))]; % center of targeting nodes
ptc  = [mean(posx([id,it])),mean(posy([id,it]))]; % center of stimulus and target
ptx(idd(1:2:end)) = 1;
pty(idd(2:2:end)) = 1;
ttx(itt(1:2:end)) = 1;
tty(itt(2:2:end)) = 1;
tx(idd(1:2:end))  = 1;
tx(itt(1:2:end))  = 1;
ty(idd(2:2:end))  = 1;
ty(itt(2:2:end))  = 1;
for i = 1:np
    dx = posx(id(i)) - pcen(1);
    dy = posy(id(i)) - pcen(2);
    prot(2*id(i)-1) = -dy;
    prot(2*id(i))   = dx;
    dx = posx(id(i)) - ptc(1);
    dy = posy(id(i)) - ptc(2);
    rot(2*id(i)-1)  = -dy;
    rot(2*id(i))    = dx;
end
for i = 1:nt
    dx = posx(it(i)) - tcen(1);
    dy = posy(it(i)) - tcen(2);
    trot(2*it(i)-1) = -dy;
    trot(2*it(i))   = dx;
    dx = posx(it(i)) - ptc(1);
    dy = posy(it(i)) - ptc(2);
    rot(2*it(i)-1)  = -dy;
    rot(2*it(i))    = dx;
end
ptx  = ptx/norm(ptx);
pty  = pty/norm(pty);
prot = prot/norm(prot);
ttx  = ttx/norm(ttx);
tty  = tty/norm(tty);
trot = trot/norm(trot);
tx   = tx/norm(tx);
ty   = ty/norm(ty);
rot  = rot/norm(rot);
for d=1:dim
    pts(idd(d:dim:end),d) = 1;
    tts(itt(d:dim:end),d) = 1;
    pts(:,d)  = pts(:,d)/norm(pts(:,d));
    tts(:,d)  = tts(:,d)/norm(tts(:,d));
end
pert = pert-pts*pts'*pert-prot*prot'*pert;
targ = targ-tts*tts'*targ-trot*trot'*targ;

% on perturbing sites
A = [pts(idd,:)';prot(idd,:)'];
naa = null(A);
Pmat = [naa';A];
dP0 = size(naa,2);
dP1 = size(A,1);
% on target sites
A = [tts(itt,:)';trot(itt,:)'];
naa = null(A);
Tmat = [naa';A];
dT0 = size(naa,2);
dT1 = size(A,1);
% on perturbing and target sites
PTmat = zeros(size(Pmat)+size(Tmat));
PTmat(1:np*dim,1:np*dim) = Pmat;
PTmat(np*dim+(1:nt*dim),np*dim+(1:nt*dim)) = Tmat;
%% Compute Mmatrix
rbb  = zeros(nnb,nnb); % boundary boundary distance
bpx  = zeros(1,nnb);
bpy  = zeros(1,nnb);
delx = zeros(1,nnb);
dely = zeros(1,nnb);
delr = zeros(1,nnb);
idx = zeros(1,4*nnb);
jdx = zeros(1,4*nnb);
val = zeros(1,4*nnb);
for n = 1:nnb
    ii = bondnb(n,1);
    jj = bondnb(n,2);
    delx(n) = posx(ii) - posx(jj);
    dely(n) = posy(ii) - posy(jj);
    if delx(n)>xbound/2   % periodic
        delx(n) = delx(n) - xbound;
    elseif delx(n)<-xbound/2
        delx(n) = delx(n) + xbound;
    end
    bpx(n) = posx(jj)+delx(n)/2;
    bpy(n) = posy(jj)+dely(n)/2;
    %if dely(n)>ybound/2  % not periodic
    %    dely(n) = dely(n) - ybound;
    %elseif dely(n)<-ybound/2
    %    dely(n) = dely(n) + ybound;
    %end
    delr(n) = sqrt(delx(n)^2+dely(n)^2);
    idx(4*(n-1)+1) = n;
    jdx(4*(n-1)+1) = 2*(ii-1)+1;
    val(4*(n-1)+1) = delx(n)/delr(n);
    idx(4*(n-1)+2) = n;
    jdx(4*(n-1)+2) = 2*(ii-1)+2;
    val(4*(n-1)+2) = dely(n)/delr(n);
    idx(4*(n-1)+3) = n;
    jdx(4*(n-1)+3) = 2*(jj-1)+1;
    val(4*(n-1)+3) = -delx(n)/delr(n);
    idx(4*(n-1)+4) = n;
    jdx(4*(n-1)+4) = 2*(jj-1)+2;
    val(4*(n-1)+4) = -dely(n)/delr(n);
end
Smatrix = sparse(idx,jdx,val,nnb,dim*num);
for n1=1:nnb
    for n2=n1:nnb
        dx = bpx(n2)-bpx(n1);
        if dx<-xbound/2
            dx = dx+xbound;
        elseif dx>xbound/2
            dx = dx-xbound;
        end
        dy = bpy(n2)-bpy(n1);
        dr = sqrt(dx^2+dy^2);
        rbb(n1,n2) = dr;
        rbb(n2,n1) = dr;
    end
end
rmin = min(min(rbb(rbb>0)));
rmax = max(max(rbb));
binn = floor((rmax-rmin)/drbin);
rbd  = rmin-drbin/2:drbin:rmin+(binn-1.5)*drbin;
if length(rbd)~=binn
    binn = length(rbd);
end
dn   = max(1,round(2.^(1:log2(nnb/2)/topnb:log2(nnb))-2.^(1-log2(nnb/2)/topnb:log2(nnb/2)/topnb:log2(nnb)-log2(nnb/2)/topnb)));
nbd  = cumsum(dn);
nbd  = nbd(2:end-1);
if topnb~=length(nbd)
    topnb = length(nbd);
end
%Tmatrix = Smatrix';
% of the weak network
bpxw = zeros(1,nnb);
bpyw = zeros(1,nnb);
delxw = zeros(1,nwk);
delyw = zeros(1,nwk);
delrw = zeros(1,nwk);
idxw = zeros(1,4*nwk);
jdxw = zeros(1,4*nwk);
valw = zeros(1,4*nwk);
for n = 1:nwk
    ii = bondwk(n,1);
    jj = bondwk(n,2);
    delxw(n) = posx(ii) - posx(jj);
    delyw(n) = posy(ii) - posy(jj);
    if delxw(n)>xbound/2
        delxw(n) = delxw(n) - xbound;
    elseif delxw(n)<-xbound/2
        delxw(n) = delxw(n) + xbound;
    end
    %if delyw(n)>ybound/2
    %    delyw(n) = delyw(n) - ybound;
    %elseif delyw(n)<-ybound/2
    %    delyw(n) = delyw(n) + ybound;
    %end
    bpxw(n) = posx(ii)+delxw(n)/2;
    bpyw(n) = posy(ii)+delyw(n)/2;
    delrw(n) = sqrt(delxw(n)^2+delyw(n)^2);
    idxw(4*(n-1)+1) = n;
    jdxw(4*(n-1)+1) = 2*(ii-1)+1;
    valw(4*(n-1)+1) = delxw(n)/delrw(n);
    idxw(4*(n-1)+2) = n;
    jdxw(4*(n-1)+2) = 2*(ii-1)+2;
    valw(4*(n-1)+2) = delyw(n)/delrw(n);
    idxw(4*(n-1)+3) = n;
    jdxw(4*(n-1)+3) = 2*(jj-1)+1;
    valw(4*(n-1)+3) = -delxw(n)/delrw(n);
    idxw(4*(n-1)+4) = n;
    jdxw(4*(n-1)+4) = 2*(jj-1)+2;
    valw(4*(n-1)+4) = -delyw(n)/delrw(n);
end
SMatW = sparse(idxw,jdxw,valw,nwk,dim*num);
MmatW = SMatW'*SMatW;

nl = ceil(nnb/lmem);
pcstt = zeros(1,nl*lmem);
    
pos = [posx;posy];
ids = {idd,nid,itt,nit,idt,nidt};
mats= {Pmat,Tmat,PTmat};


file =[];

for cc=Minit:Mend
      
        state  = configfull(cc,:);  
        
        Mmat  = kweak*MmatW+Smatrix'*diag(state)*Smatrix;

        [pcost,~] = comppcost3PB(Mmat,pert,targ,allos,pos,ids,mats);
         
        state1 = config(cc,:);
        N=length(state1);
        distseq=zeros(N,1);

        for i=1:N
         if state1(i) == consensus(i)
          distseq(i)=0;
          else
          distseq(i)=1;
          end
         end

        dist=sum(distseq(:));
        
        file =[file; dist, pcost];
end

        dirc  = './';
        file1nm = 'prova_large';
        xname = sprintf('%d', xl);
        yname = sprintf('%d', yl);
        nsname = sprintf('%d',nsp);
        ccname = sprintf('%d',cc);
        FILE1 = [dirc,file1nm '_' xname '_' yname '_' nsname '_' ccname '_PB.dat'];

        dlmwrite(saveName,file,'\t')

end
