%% Project parent code Uses RAM only
clc;clear;close all

% energy distribution maxwellian isotropic function

phyval=calphy(struct('b',0.31e-4/4.5^3,'ne',1)); %for the upper and lower bound.

Hin = 600d3; %900 km
E0 = 3d3; % eV -> kT for chi distribution
Esig = 250; % eV sigma for mu

Emaxwell = sum((1/3*normrnd(E0,Esig,3,100)).^2,1); %Chi squared dist gives maxwell distribution
hist(Emaxwell,100)
xlabel('Energy eV')
ylabel('Counts')

info = cell(1,1d5);
Etot = E0;
create = [];
% for n = 1:numel(r)
%     vo = precip1(r(n),Hin);
%     info{n} = vo;
%     
% end
E0maxwell = cell(1,length(Emaxwell));

for i = 1:numel(Emaxwell)
    i
    v0.gen0 = [Emaxwell(i);Hin];
    

    k = 0;
    
    while Etot >25
        % V =cell(length(V{1}.vo.Emat),1);
        eval(['temp = v0.gen' num2str(k) ';'])
        
        genx = [];
        tempk = [];
        [l,w] = size(temp);
        for n = 1:w
            vn = precip1(temp(1,n),temp(2,n)); % generate simulation
            genx = [genx vn.part]; % next generation of particles
        end
        
        if isempty(genx) == 1 % check if no particles were generated
            break
        end
        
        eval(['v0.gen' num2str(k+1) '=genx;']) % save particles to structure

        Echeck = genx(1,:) > 25; % check which new particles have more than 25 eV
        
        if sum(Echeck) == 0     % if all particles have less than 25 eV break the loop
            break
        end
        
        k = k+1;
    end

    E0maxwell{i} = v0; %after one primary particle, save to workspace in a cell.
end
%% extract data and plot so far
% load('save/genall.mat')

genmat = [];
enmat = [];
for m = 1:length(E0maxwell)
    temp = E0maxwell{m};
    
    for n = 1:length(temp)
%         tmpexist = eval(['exist(''temp.gen' num2str(n) ''')']);
%         if tmpexist == 0;
%         else
            eval(['genmat = [genmat temp.gen' num2str(n) '(2,:)];']);
            eval(['enmat = [enmat temp.gen' num2str(n) '(1,:)];']);
%         end
    end
    
end

[Npart,partcent] = hist(genmat,150);

figure(2)
scatter(Npart,partcent,'.')
set(gca,'xscale','log')
% xlim([0 1d5])
ylim([0 600d3])
ylabel 'Altitude [m]'
xlabel ('Ionization Rate [cm^-3 s^-1]')
title 'Altitude vs. Ionization Rate'
saveeps('ionrate3.png')

Emaxsum = sum(Emaxwell);
enmatsum = sum(enmat);

dif = enmatsum-Emaxsum;

% %%
% bin = histcounts(create,1000);
% 
% figure(1)
% scatter(1:1000,bin,'.')
% % set(gca,'view',[90 -90])
% %%
% alpha = rand(1000,1)*2*pi;
% nalpha = (1+cos(alpha).^2);
% figure(2)
% scatter(alpha,nalpha,'.')