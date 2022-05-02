%% Viscous Flow in Europa's Ice Shell
clear variables
rhow=1050;
rhoi=920;
drho=rhow-rhoi;
g=1.315;
Tm=270;
Ts=100;
r=1535*10^3;
tconv=86400*365.25;
%% mesh:
Nz=101;
Ntheta=101;
Nt=50000;

H=zeros(1,Ntheta);Hd=zeros(Nt,Ntheta);dHdtheta=zeros(1,Ntheta); dHdt=zeros(1,Ntheta);
u=zeros(Nz,Ntheta);q=zeros(1,Ntheta);
A=zeros(Nz,3);T=zeros(Nz,1);eta=zeros(Nz,Ntheta);
dz=zeros(1,Ntheta); %grid spacing

ztdim=(rhow-2*rhoi)/(2*rhow); %z=0 is neutral buoyancy
zdim=linspace(-rhoi/rhow,ztdim,Nz);%dimensionless vertical coordinate (z/H)
dzdim=abs(zdim(2)-zdim(1)); %zt is transitional depth
theta=linspace(0,pi,Ntheta); %co-latitude
ds=theta(2)-theta(1);
dt=tconv;
%% initial perturbation:
Hmin=2*10^3;
Havg=25*10^3;
Heq=0.5*(Hmin+Havg);

dH=Havg-Hmin;
dtheta=pi/2;

for i=1:Ntheta
    if theta(i)<=pi/2+dtheta && theta(i)>=pi/2-dtheta
H(i)=Hmin+dH/dtheta*abs(theta(i)-pi/2); %initial shell thickness
    else
H(i)=Hmin+dH;
    end
end
 dHdtheta(1:(Ntheta-1)/2)=-dH/dtheta; %basal slope
 dHdtheta((Ntheta-1)/2+2:Ntheta)=dH/dtheta;
%  dHdtheta(:)=dH/dtheta;
 
Hd(1,:)=H; %prescribe initial thickness
dz(:)=dzdim*Hd(1,:);
%% rheology:
eta0=10^12; %basal viscosity [diffusion creep]
%A0disCOLD=4*10^4*(10^24);%[Pa^-4/s] n=4: dislocation creep
%A0disWARM=6*10^28*(10^24);%[Pa^-4/s] n=4: dislocation creep
%A0GBSCOLD=3.9*10^(-3)*10^(6*1.8);%[Pa^-1.8/s] n=1.8: grain boundary sliding
%A0GBSWARM=3*10^26*10^(6*1.8);%[Pa^-1.8/s] n=1.8: grain boundary sliding
R=8.31; %Gas constant [J/mol/K]
A0=zeros(1,3);
n=[1, 1.8, 4];
Q=[60, 49, 60].*10^3;%activation energy [J/mol]

% creep=[1 2 3];
% [1]= newtonian diffusion creep
% [2]= grain boundary sliding, T<=259K
% [3]= dislocation creep,      T<=259K

for i=1:3
    A0(i)=1/(2*eta0)*... %effective viscosity
    (rhoi*g*dzdim*mean(H)*dH/dtheta/r*drho/rhow)^(1-n(i));
for j=1:Nz
   T(j)=Tm+(Ts-Tm)*(zdim(j)-zdim(1)); %temperature

   A(j,i)=A0(i)*exp(-Q(i)/R*(1/T(j)-1/Tm)); %creep parameter   
end
end
eta(1,:)=eta0;
creep=1
%% Simulation    

for k=2:Nt
    for j=Nz-1:-1:2
        for i=2:Ntheta-1
        
u(j-1,i)=u(j+1,i)-2*dz(i)*...%meridional ice velocity
    (2*A(j,creep)*...
    (rhoi*g*(zdim(j)-zdim(1))*Hd(k-1,i)*drho/rhow*abs(dHdtheta(i))/r)^n(creep))...
    *sign(dHdtheta(i));

q(i)=Hd(k-1,i)*u(1,i)/2; %ice flux

% eta(j,i)=1/(2*A(j,creep))*... %effective viscosity
%    (rhoi*g*(zdim(j)-zdim(1))*Hd(k-1,i)*abs(dHdtheta(i))/r*drho/rhow)^(1-n(creep));
        end
    end
for i=2:Ntheta-1 %thinning rate
     dHdt(i)= -1/r*((q(i+1)-q(i-1))/(2*ds) + ...
    q(i)*cot(theta(i))); 
end
dHdt(1)=dHdt(2); dHdt(end)=dHdt(end-1);

for i=2:Ntheta-1 % thickness for next timestep
Hd(k,i)=Hd(k-1,i)+dt*dHdt(i);
end
Hd(k,1)=Hd(k,2); Hd(k,end)=Hd(k,end-1);

for i=1:Ntheta
dz(i)=dzdim*Hd(k,i); %grid spacing for next timestep
if i>1 && i<Ntheta
        %thickness gradient for next timestep:
    dHdtheta(i)=(Hd(k,i+1)-Hd(k,i-1))/(2*ds); 
end  
end
dHdtheta(1)=dHdtheta(2); dHdtheta(end)=dHdtheta(end-1);

end

%% figures
% figure(1)
% subplot(2,1,1)
% hold on
% grid on
% plot(theta./pi,H./10^3,'b','LineWidth',2)
% plot(theta./pi,ones(size(H)).*Heq./10^3,'c--','LineWidth',2)
% xlabel('\theta [pi]')
% ylabel('Shell Thickness')
% set(gca, 'YDir', 'Reverse')
% ylim([0 Hmin+dH]./10^3)
% 
% subplot(2,1,2)
% hold on
% grid on
% plot(theta./pi,dHdt./10^3.*tconv,'LineWidth',2)
% xlabel('\theta [pi]')
% ylabel('dHdt [km/yr]')
% set(gca, 'YDir', 'Reverse')
%%

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
%%

% figure(3)
% hold on
% grid on
% plot(theta(2:Ntheta-1)./pi,q(2:Ntheta-1).*tconv)
% xlim([theta(2) theta(Ntheta-1)]./pi)
% ylabel('flow magnitude [m^2/yr]')
% xlabel('\theta [pi]')
% title('Volumetric Flow Rate')
%%
% figure(4)
% hold on
% imagesc(theta./pi,zdim-zdim(1),-u.*tconv)
% set(gca,'YDir','normal')
% ylim([0, ztdim-zdim(1)])
% xlim([theta(2) theta(Ntheta-1)]./pi)
% xlabel('\theta [pi]');
% ylabel('height in shell');
% title('Equatorial Ice Velocity')
% clrbrv=colorbar;
% clrbrv.Label.String = 'u_{\theta} [m/yr]';
%%
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
%%
num=10;
figure(7)
hold on
grid on
%plot(theta,Hd(1,:)./10^3,'b-','LineWidth',2)
for k=1:num
    plot(theta./pi,Hd(k*Nt/num,:)./10^3)
end
xlim([theta(1) theta(Ntheta)]./pi)
xlabel('\theta [pi]')
ylabel('Shell Thickness [km]')
set(gca, 'YDir','Reverse')

%%
%H12full=Hd(1:100:50000,:);