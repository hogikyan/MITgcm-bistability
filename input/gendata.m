clear all

fmt='real*8';
accuracy='real*8';
Ieee='b';


%---> Parameters

nx=105;
ny=70;
nz=29;
onedegx=nx/210;
onedegy=ny/140;

NyNorth=ny;

Onebasin = 0;

lati=2*((1:ny)-(ny+1)/2);

time=inf;

%exp=cellstr(['LGM_single_ice/output/present_idiff50']);
%dxc=rdmds(['../../',char(exp(1)),'/DXC']);
%dyc=rdmds(['../../',char(exp(1)),'/DYC']);% constant!

latc=ones(nx,1)*[-69:2:69];
dyc=6.37e6*2*pi/360*ones(nx,ny);
dxc=6.37e6*2*pi/360*ones(nx,ny).*cosd(latc);

%================================================
%=============>   TOPOGRAPHY    <================
%================================================

topo = zeros(nx,ny);
topo(:,:) = -4000. ;

%%south
%topo(:,1) = 0.;
%%north
%%topo(:,ny) = 0.;
%topo(:,NyNorth:ny) = 0.;
%%America  
%topo(1 ,23:ny) = 0.;
%topo(end ,23:ny) = 0.;
%%Antarctic peninsula  
%%topo(1 ,2:6) = 0.;
%%Scotia Ridge
%topo(1 ,2:22) = -3000.;
%topo(end ,2:22) = -3000.;

%W = 25;
%h1=1500;
%h2=3500;
%AMP=8;

%Scotia Ridge
%for j=2:18
%x=  AMP*sin((j+1)*pi/(17))-2;
%topo(:,j) = -4000+h1.*exp(-(x-(1:nx)).^2/W);
%end
%for j=2:11
%x=  AMP*sin((j+1)*pi/(17))-2;
%topo(40:end,j) = -4000+h1.*exp(-(nx+x-1-(40:nx)).^2/W);
%end
%x=  AMP*sin((11+1)*pi/(17))-2;
%topo(40:end,12) = -4000+h1.*exp(-(nx+x-1-(40:nx)).^2/W);
%topo(40:end,13) = -4000+h1.*exp(-(nx+x-1-(40:nx)).^2/W);
%x=  AMP*sin((9+1)*pi/(17))-2;
%topo(40:end,14) = -4000+h1.*exp(-(nx+x-1-(40:nx)).^2/W);

%south
topo(:,1) = 0.;
%north
topo(:,ny) = 0.;
%America
topo(1 ,8:ny) = 0.;
topo(nx,8:ny) = 0.;
%Europe+Africa
topo(30,18:ny) = 0. - 4000*Onebasin;
%Asia
topo(30:nx,ny-8:ny) = 0.;
%Scotia Ridge
topo(1 ,2:7) = -3000.;
topo(nx,2:7) = -3000.;

mask = topo;

for j=1:ny
    for i=1:nx
        if mask(i,j) == 0
            mask(i,j) = 0. ; 
        else
            mask(i,j) = 1. ; 
        end
    end
end

fid=fopen('topog.bin','w','b'); fwrite(fid,topo,'real*4'); fclose(fid);

figure(1)
clf
imagesc(topo); colorbar

%================================================
%==============>      Wind      <================
%================================================

wind  = zeros(nx,ny);
tau0=10; % notice that this is actually peak wind-speed now

for j=1:20 
	y=j-0.6;
	wind(:,j)=1.007*tau0*(sin(3.5*pi*y/(ny-1)));
end
for j=21:30 
	y=j-20.4;
	wind(:,j)=-0.65*tau0*(sin(3.5*pi*y/(ny-1)));
end
for j=31:40 
	wind(:,j)=-0.65*tau0;
end
for j=40:50 
	y=j-50.1;
	wind(:,j)=0.65*tau0*(sin(3.5*pi*y/(ny-1)));
end
for j=51:70 
	y=j-60;
	wind(:,j)=0.58*tau0*(cos(3.5*pi*(y)/(ny-1)));
end

for j=23:47 
	y=j-35.25;
	wind(:,j)=wind(:,j)+0.56*tau0*(cos(3*pi*y/(ny+0.25)))^2;
end

for i=1:nx
  wind(i,:)=smooth(wind(i,:),5);
end

wind(:,1)= 0;
wind(:,70)= 0;

fid=fopen('wind_x.bin','w','b'); fwrite(fid,wind,'real*4'); fclose(fid);

figure(2)
clf
subplot(2,1,1)
plot(wind(2,:))
hold on


% meridional (katabatic) wind:
Vmax=2; % peak meridional wind-speed (m/s)

vkata=zeros(nx,ny);
for j=2:4 % decay from vkata(:,2)=Vmax to vkata(:,4)=0
    for i=1:nx % turn on oscillate in x as sin(x)
        vkata(i,j)=Vmax*abs(sin(pi*(i-0)/nx))*(cos(pi*(j-2)/4));	  
    %vkata(:,j)=Vmax*(cos(pi*(j-2)/4));
    end
end

plot(mean(vkata(:,:),1),'r')

figure(6)
clf
imagesc(vkata)

fid=fopen('wind_y.bin','w','b'); fwrite(fid,vkata,'real*4'); fclose(fid);

%================================================
%==============>  Wind Speed    <================
%================================================

wspeed  = zeros(nx,ny);

wspeed=(wind.^2+8.^2).^0.5;

wspeed(:,1)= 0;
wspeed(:,70)= 0;

fid=fopen('wspeed.bin','w','b'); fwrite(fid,wspeed,'real*4'); fclose(fid);

figure(12)
clf
subplot(2,1,1)
plot(wspeed(2,:))
subplot(2,1,2)
c_drag1=0.0027;c_drag2=0.000142;c_drag3=0.0000764;
wsm=wspeed;wsm(wsm<0.5)=0.5;atmrho=1.2;
c_drag=c_drag1./wsm+c_drag2+c_drag3.*wsm;
plot(atmrho*c_drag(2,:).*wspeed(2,:).*wind(2,:));
hold on

%=====================================================
%===============>  Air temperature  <=================
%=====================================================


T_air = zeros(nx,ny);

theta = ny-3; %width

for i=1:nx
 for j=2:ny-1
     y=j-(ny+1)/2;
     if y<0;
         T_air(i,j) = mask(i,j)*(273.15+(37.5-3)*cos(pi*y/theta)-5.5);
     else
         T_air(i,j) = mask(i,j)*(273.15+(29.0+0)*cos(pi*y/theta)-0.0);
     end
 end
end

fid=fopen('T_air.bin','w','b'); fwrite(fid,T_air,'real*4'); fclose(fid);

figure(3)
clf
plot(T_air(10,2:end-1)-273.15);
hold off;

%====================================================================
%===================> Radiative flux <============================
%====================================================================


LWflx=0*T_air;
fid=fopen('LW_flux.bin','w','b'); fwrite(fid,LWflx,'real*4'); fclose(fid);

SWflx=0*T_air;
fid=fopen('SW_flux.bin','w','b'); fwrite(fid,SWflx,'real*4'); fclose(fid);


%=====================================================
%==================>    E - P    <====================
%=====================================================
empconst=0.5/3.1536e7; 

theta1 = ny-3;

emp = zeros(nx,ny);

for i=1:nx % Southern Ocean
 for j=2:19
     y=j-(ny+1)/2;
     emp(i,j) = -empconst*mask(i,j)*(cos(1.9*pi*abs(2*(y-3)/theta1)^0.75))*(2.6*cos(0.8*pi*abs((y-3)/theta1)^1.0))+0.2*10^-8;
 end
 for j=16:25
     y=j-(ny+1)/2;
     emp(i,j) = -empconst*mask(i,j)*(cos(1.9*pi*abs(2*(y-2)/theta1)^0.75))*(2.5*cos(0.8*pi*abs((y-3)/theta1)^0.9))+0.2*10^-8;
 end
end

for i=30:nx % Indo-Pacific
 for j=26:ny-1
     y=j-(ny+1)/2;
     emp(i,j) = -empconst*mask(i,j)*(cos(2*pi*abs(2*(y-3)/theta1)^0.5))*(2.5*cos(1*pi*abs(y/theta1)^0.8))-0.8*10^-8;
    if y<0
       emp(i,j) = -empconst*mask(i,j)*(cos(1.9*pi*abs(2*(y-2)/theta1)^0.75))*(2.2*cos(0.8*pi*abs((y-3)/theta1)^1.5))+0.2*10^-8;
       if y>-3
          emp(i,j) = -empconst*mask(i,j)*(sin(1.9*pi*abs(2*(y-2)/theta1)^0.75))*(0.5*cos(0.8*pi*abs((y-3)/theta1)^1.5))+0.2*10^-8;
       end
    end
 end
end

for i=1:30 % Atlantic
 for j=26:ny-1
     y=j-(ny+1)/2;
     emp(i,j) = -empconst*mask(i,j)*(cos(2.2*pi*abs(2*(y-2)/theta1)^0.6))*(2*cos(1.0*pi*abs(y/theta1)^1));
    if y<0
       emp(i,j) = -empconst*mask(i,j)*(cos(1.5*pi*abs(2*(y-5)/theta1)^0.65))*(2.2*cos(0.8*pi*abs((y-3)/theta1)^1.5))+0.2*10^-8;
    end
 end
end

% Add diagnosed flux from surface salinity relaxation
addpath /project2/mfj/hogikyan/MITgcm_01_20_2026/utils/matlab
S0=34.8;
diag=rdmds('InputsurfDiag',0047085000);
sflux=diag(:,:,3)/S0/10^3; % to convert to the same unit as emp
emp=emp+sflux;

emp=emp-mean(mean(emp.*dxc.*mask))/mean(mean(dxc.*mask));
emp=emp.*mask;

% CC-scaling amplification
emp=emp*((1.07)^2);

%=====================================================
%===============>    Evaporation    <=================
%=====================================================

evap=0.0*emp;
evap(emp>0)=emp(emp>0);

fid=fopen('evap.bin','w','b'); fwrite(fid,evap,'real*4'); fclose(fid);

%=====================================================
%===============>   Precipitation   <=================
%=====================================================

precip=0.0*emp;
precip(emp<0)=-emp(emp<0);

snowprecip=0*precip;
snowprecip(T_air<273.15)=precip(T_air<273.15);

fid=fopen('precip.bin','w','b'); fwrite(fid,precip,'real*4'); fclose(fid);
fid=fopen('snowprecip.bin','w','b'); fwrite(fid,snowprecip,'real*4'); fclose(fid);




figure(7)
clf
subplot(3,1,1)
plot(lati([2,end-1]),[0,0],'--k');hold on
plot(lati(2:NyNorth-1),emp(10,2:NyNorth-1));hold on;
subplot(3,1,2)
plot(lati(2:NyNorth-1),precip(10,2:NyNorth-1));hold on;
plot(lati(2:NyNorth-1),snowprecip(10,2:NyNorth-1));
subplot(3,1,3)
plot(lati(2:NyNorth-1),evap(10,2:NyNorth-1));

%=====================================================
%============>  Sea surface salinity  <===============
%=====================================================

theta1 = ny-3;
S0=34.8;

SSS = zeros(nx,ny);


for i=1:nx %Southern Ocean and Atlantic
    for j=2:ny-1
        y=j-(ny+1)/2;
        SSS(i,j) = -1.5*mask(i,j)*(cos(2*pi*(abs(2*(y-2)/theta1)^0.5)))+0.7;
        if y>=23;
           SSS(i,j) = 0.1;
        end
    end
end
for i=31:nx %Indo-Pacific
    for j=21:49
        y=j-(ny+1)/2;
        SSS(i,j) = 1.*mask(i,j)*(cos(2*pi*(abs(2*(y-19)/theta1)^1.0)))+0.25;
        SSS(i,j) = SSS(i,j) + 1*(0.5*y-abs(y))/ny;
    end
    for j=49:ny-1
        y=j-(ny+1)/2;
        SSS(i,j) = -2.25*mask(i,j)*(cos(2*pi*(abs(2*(y-2)/theta1)^0.5)))-1.25;
    end
end
SSS=SSS+S0;
SSS=SSS.*mask;

fid=fopen('SSS.bin','w','b'); fwrite(fid,SSS,'real*4'); fclose(fid);

figure(8)
clf
plot(lati(2:NyNorth-1),SSS(10,2:NyNorth-1));hold on;
plot(lati(2:NyNorth-1),SSS(40,2:NyNorth-1));
ylim([S0-3,S0+3]);xlim([-70,70]);hold off;
