%% Project parent code uses RAM and physical memory
clc;clear;close all

% energy distribution maxwellian isotropic function

phyval=calphy(struct('b',0.31e-4/4.5^3,'ne',1)); %for the upper and lower bound.

Hin = 900d3; %900 km
E0 = 2d3; % eV -> kT for chi distribution
Esig = 1000; % eV sigma for mu

Emaxwell = sum((normrnd(E0,Esig,3,10000)).^2,1); %Chi squared dist gives maxwell distribution
hist(Emaxwell,100)

info = cell(1,1d5);
Etot = E0;
create = [];
% for n = 1:numel(r)
%     vo = precip1(r(n),Hin);
%     info{n} = vo;
%     
% end

alpha = rand(1000,1)*2*pi;
nalpha = (1+cos(alpha).^2)*pi;
figure(2)
histogram(nalpha,100)

for i = 1:numel(Emaxwell)
    vo.genx = [Emaxwell(i);Hin];
    savename = 'save/gen0.mat';
    save(savename,'vo');

    k = 0;

    while Etot >25

        loadname = ['save/gen' num2str(k) '.mat'];
        vi = load(loadname);

        % V =cell(length(V{1}.vo.Emat),1);
        genx = [];
        [l,w] = size(vi.vo.genx);
        for n = 1:w
            vn = precip1(vi.vo.genx(1,n),vi.vo.genx(2,n));
            genx = [genx vn.part];
            temp.genx = genx;
        end

        if isempty(genx) == 1
            break
        end

        vo.genx = temp.genx;
        savename = ['save/gen' num2str(k+1)];
        save(savename,'vo');


        Echeck = vo.genx(1,:) > 25;

        if sum(Echeck) == 0
            break
        end

        k = k+1;
    end



    for m = 1:k
        loadname = ['save/gen' num2str(m) '.mat'];
        vf = load(loadname);
        create = [create vf.vo.genx(2,:)];
    end
end

%%
bin = histcounts(create,1000);
figure(1)
histogram(create)
set(gca,'view',[90 -90])



alpha = rand(1000,1)*2*pi;
nalpha = (1+cos(alpha).^2);
figure(2)
scatter(alpha,nalpha,'.')