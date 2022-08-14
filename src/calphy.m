function val=calphy(strin,redefine)
% References: 
% 1) F. F. Chen, Chinese Edition, p206-207
% 2) P. Gary, 1993, p5
% specify the background plasma
% then output variable parameter
% e: electron, h : proton
% w: omega 

% 9/17/2012 add dynamics of particles in a dipole field
% 12/25/2014 allow redefining constants
%   redefine variables are passed through redefine={constant1name,constant1value,constant2name,constant2value,...}
% 12/25/2014 allow changing default constants. phyval_global_
% 
hmas=1.67252e-27;
ev  =1.60210e-19;
emas=9.10908e-31;
%e0  =8.85442e-12; 
e0   =8.854187817e-12; % permittivity of vaccum C^2/J m
%ee0 =e0;
mu0 =4.0*pi*1e-7; % permeability of vacuum T^2 m^-3/J
% clight=2.99792501e8;
clight=2.997924580105029e8;  % speed of light m/s
kb  =1.38054e-23;
RE=6370d3; % Earth's radius
B0dip=0.31d-4;

if (nargin ==0) 
    strin=struct('ne',1,'b',1);
end

if (strcmp(strin,'constant'))
	  o2.hmas=hmas; 
      o2.emas=emas;
      o2.ev=ev;
      o2.e0=e0;
      o2.mu0=mu0;
      o2.clight=clight;
      o2.kb=kb;
      o2.RE=RE;
      o2.B0dip=B0dip;
      val=o2;
      return
end




global phyval_global_
if (strcmp(strin,'reset'))
    phyval_global_=[];
end


if isstruct(phyval_global_)
    % using the global variable as default
    % this can be turned off by setting phyval_global_=[]
    snames=fieldnames(phyval_global_);
    for isn=1:length(snames)
        cmd=[snames{isn},'=getfield(phyval_global_,snames{isn})'];
        eval(cmd);
    end
end






if nargin>=2 % redefine constants here
    % strin.redefine={'hmas',hmas*2}
    nchange=length(redefine)/2;
    for ich=1:nchange
        cmd=[redefine{(ich-1)*2+1},'=redefine{ich*2}'];
        disp(cmd);
        eval(cmd);
    end
end
%ee0 =e0;



rad =pi/180.;
epsilon=emas/hmas;
memp=epsilon;

val.mu0=mu0;
val.e0=e0;
val.kb=kb;
val.B0dip=B0dip;
val.RE=RE;
val.Erest_e=emas*clight^2/ev;
val.Erest_p=hmas*clight^2/ev;

% clj-RN(I)-pp7
% conversion f[s3m-6] = fac* J[cm-2 s-1 ster-1 kev-1]/(E[keV])^2
val.Hfac=5.4490e-19;
val.Hefac=8.7185e-18;
val.Ofac=1.3950e-16;
val.efac=1.6163e-25;

val.gainfac=20*log10(exp(1)); % gain in dB = gfac*val.gainfac
val.dbfac=20*log10(exp(1))*6400e3; % factor converting ki to dB/RE





if (~ isfield(strin,'ne'))
	  strin.ne=1; % m-3
end

if (~ isfield(strin,'b'))
	  if(isfield(strin,'l'))
		strin.b=0.31e-4./strin.l.^3;
	  else
	  	strin.b=1; % T
      end
end

if (~ isfield(strin,'te'))
	  strin.te=1; % ev
end

if (~ isfield(strin,'th'))
	  strin.th=1; % ev
end

ne=strin.ne;
b=strin.b;
if strin.ne<0 % wpe 
    ne=abs(strin.ne.^2)*e0*emas/ev/ev;
end
if strin.b<0 % wce
    b=abs(strin.b)*emas/ev;
end


% default ion compisition
if (~ isfield(strin,'etahe'))
    val.etah=1;
    val.etahe=0;
    val.etao=0;
else
    val.etah=strin.etah;
    val.etahe=strin.etahe;
    val.etao=abs(1-strin.etah-strin.etahe);
end

val.chfr=chfr_(val.etah,val.etahe);

if (isfield(strin,'plsmat')) % hotray distribution function
% plsmat variable is consistent with plsmat@hpdr
% temperature is specified in keV
% 9/5/2012
    hdf=strin.fv;
    indc=find(hdf(:,6)>0);
    hdf=hdf(indc,:);
    sz=size(hdf);
%     hdf=[...
% 1.0E0 0.0   T_Oc T_Oc 0.5  N_Oc 16.  1; ... %1
% 1.0E0 0.0   T_ec T_ec 0.5  N_ec 0.  -1; ... %2
% 1.0E0 0.0   T_Hec T_Hec 0.5  N_Hec 4.  1; ... %3
% 1.0E0 0.0   T_Hc T_Hc 0.5  N_Hc 1.  1; ... %4
% 1.0E0 0.0   Tp_Hh Tz_Hh 0.5  N_Hh  1. 1; ... %5
% 1.0E0 0.0   Tp_eh Tz_eh 0.5  N_eh 0. -1; ... %6
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %7
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %8
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %9
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %10
% ]; 
    
    
    val.ncomp=sz(1);
    indel=find(hdf(:,7)==0 & hdf(:,8)==-1);
    indh=find(hdf(:,7)==1 & hdf(:,8)==1);
    indhe=find(hdf(:,7)==4 & hdf(:,8)==1);
    indo=find(hdf(:,7)==16 & hdf(:,8)==1);
    val.ne=sum(hdf(indel,6));
    ne=val.ne;
    if(~isempty(indh)), val.etah=sum(hdf(indh,6))./val.ne;,end
    if(~isempty(indhe)), val.etahe=sum(hdf(indhe,6))./val.ne;,end
    if(~isempty(indo)), val.etao=sum(hdf(indo,6))./val.ne;,end
    
    for ic=1:val.ncomp
        
        val.df(ic).mas=hdf(ic,7)*hmas;
        if(hdf(ic,7)==0), val.df(ic).mas=emas;, end
        val.df(ic).q=hdf(ic,8);
        val.df(ic).n=hdf(ic,6);
        val.df(ic).Tz=hdf(ic,4); % eV
        val.df(ic).Tp=hdf(ic,3); % eV
        
        val.df(ic).wp=sqrt(val.df(ic).q.^2*ev*ev/val.df(ic).mas/e0*val.df(ic).n);
        val.df(ic).wc=val.df(ic).q*ev/val.df(ic).mas*b;
        val.df(ic).az=sqrt(2*val.df(ic).Tz*ev/val.df(ic).mas);
        val.df(ic).ap=sqrt(2*val.df(ic).Tp*ev/val.df(ic).mas);
    end
end

if (isfield(strin,'fv')) % hotray distribution function
% this fv field was used before 9/15/2012
% now will use plsmat field, where temperature is specified in keV
% rather than eV as fv field

    hdf=strin.fv;
    indc=find(hdf(:,6)>0);
    hdf=hdf(indc,:);
    sz=size(hdf);
%     hdf=[...
% 1.0E0 0.0   T_Oc T_Oc 0.5  N_Oc 16.  1; ... %1
% 1.0E0 0.0   T_ec T_ec 0.5  N_ec 0.  -1; ... %2
% 1.0E0 0.0   T_Hec T_Hec 0.5  N_Hec 4.  1; ... %3
% 1.0E0 0.0   T_Hc T_Hc 0.5  N_Hc 1.  1; ... %4
% 1.0E0 0.0   Tp_Hh Tz_Hh 0.5  N_Hh  1. 1; ... %5
% 1.0E0 0.0   Tp_eh Tz_eh 0.5  N_eh 0. -1; ... %6
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %7
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %8
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %9
% 1.0E0 0.0   1e-3 1e-3 0.5  0e6  0. -1; ... %10
% ]; 
    
    
    val.ncomp=sz(1);
    indel=find(hdf(:,7)==0 & hdf(:,8)==-1);
    indh=find(hdf(:,7)==1 & hdf(:,8)==1);
    indhe=find(hdf(:,7)==4 & hdf(:,8)==1);
    indo=find(hdf(:,7)==16 & hdf(:,8)==1);
    val.ne=sum(hdf(indel,6));
    ne=val.ne;
    if(~isempty(indh)), val.etah=sum(hdf(indh,6))./val.ne;,end
    if(~isempty(indhe)), val.etahe=sum(hdf(indhe,6))./val.ne;,end
    if(~isempty(indo)), val.etao=sum(hdf(indo,6))./val.ne;,end
    
    for ic=1:val.ncomp
        
        val.df(ic).mas=hdf(ic,7)*hmas;
        if(hdf(ic,7)==0), val.df(ic).mas=emas;, end
        val.df(ic).q=hdf(ic,8);
        val.df(ic).n=hdf(ic,6);
        val.df(ic).Tz=hdf(ic,4); % eV
        val.df(ic).Tp=hdf(ic,3); % eV
        
        val.df(ic).wp=sqrt(val.df(ic).q.^2*ev*ev/val.df(ic).mas/e0*val.df(ic).n);
        val.df(ic).wc=val.df(ic).q*ev/val.df(ic).mas*b;
        val.df(ic).az=sqrt(2*val.df(ic).Tz*ev/val.df(ic).mas);
        val.df(ic).ap=sqrt(2*val.df(ic).Tp*ev/val.df(ic).mas);
    end
end


% if isfield(strin,'dipolemotion')
%     % dynamics in dipole field
%     % in relativistic formula
%     % 
% end

nh=ne;
te=strin.te;
th=strin.th;

val.emas=emas;
val.hmas=hmas;
val.memp=emas/hmas;
val.memh=emas/hmas;
val.mhme=hmas/emas;
val.mpme=hmas/emas;
val.ev=ev;
val.c=clight;
val.clight=clight;
val.kb=kb;

val.ne=ne;
val.b=b;
val.te=te;
val.th=th;

val.wce=ev/emas*b;
val.wch=ev/hmas*b;
val.wpe=sqrt(ev*ev/emas/e0*ne);
val.wph=sqrt(ev*ev/hmas/e0*nh);
val.fce=val.wce/(pi*2);
val.fch=val.wch/(pi*2);
val.fpe=val.wpe/(pi*2);
val.fph=val.wph/(pi*2);
val.fpara=val.fpe./val.fce;
val.rho=nh*hmas+ne*emas;
val.va=b./sqrt(mu0*val.rho); % F. F. Chen p207


% lower hybrid resonance frequency here with density correction
% see lowerhybridresonance_pic.pdf
val.wlhr=(1./(val.wce.*val.wch)+1./(val.wph.^2)).^(-0.5);
val.flhr=val.wlhr/(pi*2);




val.ve=sqrt(2*te*ev/emas); % definition of thermal speed
val.vh=sqrt(2*th*ev/hmas);
val.rhoe=val.ve./val.wce;
val.rhoh=val.vh./val.wch;

% Maxwell Speed Distribution, MSD ~ v^2 exp(-v^2/a^2) where v from 0 to Inf
% http://en.wikipedia.org/wiki/Maxwell_speed_distribution
% There are three candidates for the averaged value of the MSD
% Here a is the speed of maximum MSD (Most Probable speed, v_mp)
% Mean speed <v> = sqrt(4/pi) * a
% root mean squared speed v_rms = sqrt(3/2) * a
% v_mp < <v> < v_rms
% Note that in Gary 1993, v=a/sqrt(2) wildly used
val.ve_g93=val.ve/sqrt(2);
val.vh_g93=val.vh/sqrt(2);
val.ke_g93=val.wpe./val.ve_g93;
val.kh_g93=val.wph./val.vh_g93;



% Debye length, Debye wave number
val.De=sqrt(e0*(te*ev)./ne/ev^2); % F. F. Chen p207
val.De=val.ve./sqrt(2)./val.wpe;
val.ke=1./val.De;

val.Dh=sqrt(e0*(th*ev)./nh/ev^2); % ~ (T/n) ^(1/2) q^-1 
val.Dh=val.vh./sqrt(2)./val.wph;
val.kh=1./val.Dh;

% plasma beta, plasma pressure / magnetici pressure
val.bet=(ne*te*ev+nh*th*ev)./(b.^2/2/mu0);  % F. F. Chen p207
val.bete=(ne*te*ev)./(b.^2/2/mu0); % Gary 1993, p5
val.beth=(nh*th*ev)./(b.^2/2/mu0);

% inertial length of particle species, Gary 1993, p5
val.ile=clight./val.wpe;
val.ilh=clight./val.wph;

% plasma parameter g
% S. Ichimura 1973, Basic Principles of Plasma Physics, A statistical approach
% equ (1.22) p11
val.g=1./(val.ne.*val.De.^3); 

% plasma parameter, Gamma (equ 1.19, p15)
% The Physics of Plasmas, Richard Fitzpatrick, UTexas, Austin
% farside.ph.utexas.edu/teaching/plasma/
val.Gamma=4*pi*val.ne.*val.De.^3;



function chfrout=chfr_(etaH,etaHe)
% chfr normalized to O+ gyro
% etaH=etaarr(1);
% etaHe=etaarr(2);

fc1=16; fc2=4; fc3=1; 
xx=etaHe/etaH; yy=(1-etaH-etaHe)/etaH;
% calculate cross-over frequency 
A=1+xx+yy;
B=-(fc2^2+fc3^2+xx*(fc1^2+fc3^2)+yy*(fc2^2+fc1^2));
C=fc2^2*fc3^2+xx*fc1^2*fc3^2+yy*fc1^2*fc2^2;
fcr1=sqrt((-B+sqrt(B^2-4*A*C))/2/A);
fcr2=sqrt((-B-sqrt(B^2-4*A*C))/2/A);

% calculate lcut off frequency 
A=1+xx+yy;
B=-(fc2+fc3+xx*(fc1+fc3)+yy*(fc2+fc1));
C=fc2*fc3+xx*fc1*fc3+yy*fc1*fc2;
flcut1=(-B+sqrt(B^2-4*A*C))/2/A;
flcut2=(-B-sqrt(B^2-4*A*C))/2/A;

% calculate bi-ion frequency 
A=fc1+xx*fc2+yy*fc3;
B=-(fc1*(fc2^2+fc3^2)+xx*fc2*(fc1^2+fc3^2)+yy*fc3*(fc2^2+fc1^2));
C=fc1*(fc2^2*fc3^2)+xx*fc2*(fc1^2*fc3^2)+yy*fc3*(fc2^2*fc1^2);

fbion1=sqrt((-B+sqrt(B^2-4*A*C))/2/A);
fbion2=sqrt((-B-sqrt(B^2-4*A*C))/2/A);

chfrout.fcr1=fcr1;
chfrout.fcr2=fcr2;
chfrout.flcut1=flcut1;
chfrout.flcut2=flcut2;
chfrout.fbion1=fbion1;
chfrout.fbion2=fbion2;


