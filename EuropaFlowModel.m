%% Viscous Creep in a Spherical Shell
% Jakob Kintzele, Princeton University Geosciences 
% Last Update: May 5, 2022

clear variables
%% =========== Constants =========== 
rhoi=920; %shell density
rhow=1050; %ocean density [if relevant]
drho=rhow-rhoi; %density contrast
g=1.315; %gravity
Tm=270; %shell melting temperature
Ts=100; %shell surface temperature
r=1535*10^3; %approximate shell radius 
tconv=86400*365.25; %s->yr
%% =========== Mesh =========== 

%forward-difference time scheme:
Nt=100000;%50000; %timesteps

Nz=101; %z-elements
Ntheta=51; %theta-elements
H=zeros(1,Ntheta);Hd=zeros(2,Ntheta);dHdtheta=zeros(1,Ntheta); dHdt=zeros(1,Ntheta);
u=zeros(Nz,Ntheta);q=zeros(1,Ntheta);
A=zeros(Nz,3);T=zeros(Nz,1);eta=zeros(Nz,Ntheta);
dz=zeros(1,Ntheta); 

trans_rel=0.25; % relative percentage of shell which deforms viscously
ztdim=(rhow-2*rhoi)/(rhow)*trans_rel; %z=0 is neutral buoyancy, zt is visc-elastic transition depth
zdim=linspace(-rhoi/rhow,ztdim,Nz);%dimensionless vertical coordinate (z/H)
dzdim=abs(zdim(2)-zdim(1)); %dimensionless z-spacing
theta=linspace(0,pi,Ntheta); %co-latitude
ds=theta(2)-theta(1); %theta spacing
dt=tconv/100; %timestep [tconv/100 (n=4,Ntheta=51)-tconv*5 (n=1-1.8,Ntheta=51)]
%% =========== initial thickness perturbation =========== 
Hmin=2*10^3; % minimum thickness
Havg=25*10^3; % initial polar/maximum thickness
Heq=0.5*(Hmin+Havg); % estimated equilibrium thickness

dH=Havg-Hmin; % max-min thickness 
dtheta=pi/2; % half the range of thinned region

for i=1:Ntheta %initial shell thickness
    if theta(i)<=pi/2+dtheta && theta(i)>=pi/2-dtheta %thinned region
Hd(1,i)=dH*cos(theta(i)*2)./2+Hmin+dH/2;
    else
Hd(1,i)=Hmin+dH;
    end
    dHdtheta(i)=-dH/2*sin(2*theta(i)); %thickness gradient
end

dz(:)=dzdim.*Hd(1,:); %initial grid spacing
%% =========== Rheology =========== 

% creep=[1 2 3];
% [1]= newtonian diffusion creep
% [2]= grain boundary sliding
% [3]= dislocation creep
creep=3

eta0=10^12; %basal viscosity [diffusion creep]
R=8.31; %Gas constant [J/mol/K]
A0=zeros(1,3);
n=[1, 1.8, 4];
Q=[60, 49, 60].*10^3;%activation energy [J/mol]

for i=1:3
    A0(i)=1/(2*eta0)*... %effective viscosity
    (rhoi*g*dzdim*mean(H)*dH/dtheta/r*drho/rhow)^(1-n(i));
for j=1:Nz
   T(j)=Tm+(Ts-Tm)*(zdim(j)-zdim(1)); %temperature
   A(j,i)=A0(i)*exp(-Q(i)/R*(1/T(j)-1/Tm)); %creep parameter   
end
end
eta(1,:)=eta0;

%% =========== Simulation ===========     

figure
set(gca,'YDir','Reverse')
pic=Nt/10; % Nt/ # of snapshots 

for k=1:Nt % time loop
    for j=Nz-1:-1:2 % z-loop
%------meridional velocity-------------------------------------------------       
u(j-1,:)=u(j+1,:)-2*dz.*...%
    (2*A(j,creep)*...
    (rhoi*g*(zdim(j)-zdim(1)).*Hd(1,:).*drho/rhow.*abs(dHdtheta)/r).^n(creep))...
    .*sign(dHdtheta);
%boundary conditions [if necessary]:
u(:,1)=0; u(:,end)=0; u(:,(Nz-1)/2+1)=0; 
%areal flux: 
q=Hd(1,:).*trans_rel.*u(1,:);   
%------effective viscosity------------------------------------------------- 
% eta(j,i)=1/(2*A(j,creep))*... 
%    (rhoi*g*(zdim(j)-zdim(1))*Hd(k-1,i)*abs(dHdtheta(i))/r*drho/rhow)^(1-n(creep));

    end
%-----thinning rate--------------------------------------------------------
for i=2:Ntheta-1 
     dHdt(i)= -1/r*((q(i+1)-q(i-1))/(2*ds) + q(i)*cot(theta(i))); 
end
dHdt(1)=dHdt(2); dHdt(end)=dHdt(end-1);                                     %Hd(k+1,i)=Hd(k-1,i)+2*dt*dHdt(i); 
                                                                            %Hd(k,1)=Hd(k,2); Hd(k,end)=Hd(k,end-1);                                                                            %dz(i)=dzdim*Hd(k,i); %grid spacing for next timestep
%-----timetepping----------------------------------------------------------
Hd(2,:)=Hd(1,:)+dt.*dHdt;
dz=dzdim.*Hd(2,:); 

%-----break loop for NaN values--------------------------------------------
    if any(isnan(Hd(2,:)), 'all')
    sprintf('NaN H Values')
    break
    end

%-----thickness gradient---------------------------------------------------
for i=2:Ntheta-1 
dHdtheta(i)=(Hd(2,i+1)-Hd(2,i-1))/(2*ds);                                    %dHdtheta(i)=(Hd(k,i+1)-Hd(k,i-1))/(2*ds); 
end  
dHdtheta(1)=dHdtheta(2); dHdtheta(end)=dHdtheta(end-1);

%-----plot snapshots-------------------------------------------------------
    if mod(k,pic)==0
hold on
plot(theta./pi,Hd(2,:)./10^3)
drawnow
    end
%-----reset thickness vector-----------------------------------------------
Hd(1,:)=Hd(2,:); 
end

%% =========== figures ============

%% ===Rheology===
% figure(2)
% 
% subplot(1,2,1)
% grid on
% hold on
% p1=plot(log10(A(:,1)/A0(1)),zdim-zdim(1),'b-','LineWidth',2);
% p2=plot(log10(A(:,2)/A0(2)),zdim-zdim(1),'r-','LineWidth',2);
% %p3=plot(log10(A(:,3)/A0(3)),zdim,'g-');
% legend([p1 p2], {'Diffusion','Dislocation, GBS'},'Location','SW')
% xlabel('log_{10}(A/A0)')
% ylabel('dimensionless shell height')
% ylim([0 zdim(end)-zdim(1)])
% title('Creep Factor')
% 
% subplot(1,2,2)
% grid on
% hold on
% plot(T,zdim-zdim(1),'k-','LineWidth',2)
% xlabel('T [K]')
% ylim([0 zdim(end)-zdim(1)])
% title('Temperature')

% figure(5)
% imagesc(theta./pi,zdim-zdim(1),log10(eta))%,[log10(eta0) log10(1000*eta0)])
% set(gca,'YDir','normal')
% ylim([0, ztdim-zdim(1)])
% xlim([theta(1) theta(Ntheta)]./pi)
% xlabel('effective viscosity [Pa s]')
% xlabel('\theta [pi]');
% ylabel('height in shell');
% title('Effective Ice Viscosity')
% clrbr=colorbar;
% clrbr.Label.String = 'log_{10}(\eta)';
%% ===velocity field===
% figure(4)
% hold on
% imagesc(theta./pi,zdim-zdim(1),u.*tconv)
% set(gca,'YDir','normal')
% ylim([0, ztdim-zdim(1)])
% xlim([theta(2) theta(Ntheta-1)]./pi)
% xlabel('\theta [pi]');
% ylabel('height in shell');
% title('Ice Velocity')
% clrbrv=colorbar;
% clrbrv.Label.String = 'u_{\theta} [m/yr]';

%% ===Thickness===
% num=10;
% figure(1)
% hold on
% grid on
% plot(theta./pi,H./10^3,'b-','LineWidth',2)
% %for k=1:num
%     plot(theta./pi,Hd(3,:)./10^3)
% %end
% xlim([theta(1) theta(Ntheta)]./pi)
% xlabel('\theta [pi]')
% ylabel('Shell Thickness [km]')
% set(gca, 'YDir','Reverse')
