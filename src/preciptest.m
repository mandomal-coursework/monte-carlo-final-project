clc
clear;close all

hin = 900d3; %starting height
Ein = 2d3;


npart = zeros(1,1000);
nsteps = zeros(1,1000);
Epart = [];

for m = 1:1000

%% ionization calculation
% GE = 1;
% TE = 1;
% 
% Eion = 1d3;
% 
% Esmax = (Ep-Eion)/2;
% 
% rmax = atan((Esmax-TE)/GE);

% Initialize the variables for the loop


En = Ein;
hn = hin;
% 
Emat = [];
hmat = [];
part = [];
k = 1; %while loop counter
n = 1;
while En > 0.25
    Emat(k) = En;         % save the energy
    hmat(k) = hn;         % save the height
    
    if En <0.25
        break
    elseif hn < 0
        break
    end
    
    choice = rand;
    
    ionp = (0.8)*exp(-hn/100d3);
    elap = 1-ionp;
    
    if choice <= ionp % check if ionized
        dE = rand*25;% random energy change
        part(1,n) = dE;   % this is the energy of an ionized particle
        part(2,n) = hn;   % this is the height that electron is created at
        n = n+1;          % count up to creat artificial loop of generated particles
 
    elseif choice>= elap   % if not ionized check if excited
        dE = 0;
    else
        dE = rand*20;     % random energy change
    end
    
    
    if dE > En
        break
    end
    
    if isempty(part);
    else
        sumcheck = sum(part(1,:));
        if sumcheck > Ein
            break
        end
    end
    
    
    dte = -log(rand)*10; % change the height it falls at again
    hn = hn-dte;
    
    En = En - dE;         % change the energy by how much the new ionized particle has
    

    
    k = k+1; % increase the counter to keep track of the particle
end

% output housekeeping
o.Emat = Emat;
o.hmat = hmat;
o.k = k;
o.part = part;

[l,w] = size(part);

if isempty(part)
else
    Epart = [Epart part(1,:)];
end
npart(m) = w;
nsteps(m) = k;
end
sum(npart)

if isempty(part)
else
    sum(part(1,:))
end
    % sum(nsteps)/1000

