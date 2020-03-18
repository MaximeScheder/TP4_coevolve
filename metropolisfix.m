function [grid,gridint,dEtot,M,E] = metropolisfix(N,kT,J,h,t,grid)

% Pick a sequence of random spins (with a linear index)
spin = randi(N,t,1);

dEtot=0;
E=0;
M=0;

% Evolve the system for a fixed number of steps t
for i=1:(t-1),
 for j=2:t,  
    s = spin(i);  % location of ith spin
    r = spin(j);  % location of jth spin
  
    gridp= 1-grid(:); % NOTE: here I need to invert the configuration because the ACE algorithm runs well on configurations where 0 amd 1 have been exchanged thus J and  h refer to this situation  
 
   if abs(r-s)>0,
    if (grid(r)==0 && grid(s)==1) || (grid(s)==0 && grid(r)==1),

    gridm=grid;
    gridpp= 1-grid(:);
    gridpp(s) = gridp(r);
    gridpp(r) = gridp(s);
        
    energyp = isingenergy(gridp,J,h);
    energypp = isingenergy(gridpp,J,h);

    dE = energyp - energypp; %  Calculate the change in energy of flipping the s-spin. Note that the 'flip' is (1-gridp(s)) for 0/1 sequences
    
    % Calculate the transition probability
    p = exp(dE/kT);

   % Decide if a transition will occur depending on the value of p
   if rand <= p, 
   grid(s) = gridm(r);
   grid(r) = gridm(s); 
   end
 else dE=0;
 end
 else dE=0;
 end
  if i==(t-2)/2,
   if j==t/2,
     gridint=grid;
   end
  end

 dEtot=[dEtot;dE];
 M = [M;sum(grid)/numel(grid)];
 E = [E;isingenergy(gridp,J,h)];

 end

 end

end
