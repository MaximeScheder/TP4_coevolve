function E = non_linear_energy(grid, h, V)
% Funciton that return the non linear energy of the model
% parameters :  grid - sequence of 0/1 representing the protein
%               h - inferred field
%               V - non linear function V(E)

h = reshape(h, [], 1);
E = -(transpose(h)*grid);
E = V(E);
end