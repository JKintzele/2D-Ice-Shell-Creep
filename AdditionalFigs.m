%% Additional Figures

theta13=linspace(0, pi, 101);
%%
num=4;
figure(7)
subplot(3,1,1)
hold on
grid on
plot([0 1], [4 4], 'c--','LineWidth',1)
plot(theta./pi,H12b(1,:)./10^3,'k--','LineWidth',2)
for k=1:num
    plot(theta./pi,H12b(k*80/num,:)./10^3)
end
set(gca, 'YDir','Reverse')

subplot(3,1,2)
hold on
grid on
plot([0 1], [4 4], 'c--','LineWidth',1)
plot(theta./pi,H12b(1,:)./10^3,'k--','LineWidth',2)
for k=1:num
    plot(theta13./pi,Hd13(k*250000/num,:)./10^3)
end

xlim([theta(1) theta(Ntheta)]./pi)
ylim([Hmin Havg]./10^3)

ylabel('Shell Thickness [km]')
set(gca, 'YDir','Reverse')

subplot(3,1,3)
hold on
grid on
plot([0 1], [4 4], 'c--','LineWidth',1)
plot(theta./pi,H12b(1,:)./10^3,'k--','LineWidth',2)
for k=1:num
    plot(theta13./pi,H14i(k*500/num,:)./10^3)
end
xlabel('\theta [pi]')
set(gca, 'YDir','Reverse')

%%
t12=linspace(0,50000,500);
t13=linspace(0,250000,500);
%%
figure(8)
subplot(1,2,1)
hold on 
grid on
plot([0 250], [4 4], 'c--','LineWidth',1)
%p1=plot(t12./10^3,Hd12(1:100:end,251)./10^3,'b-','LineWidth',2);
plot(t13./10^3,H12b(:,151)./10^3,'b-','LineWidth',2);
plot(t13./10^3,Hd13(1:500:end,51)./10^3,'r-','LineWidth',2);
plot(t13./10^3,H14i(:,51)./10^3,'g-','LineWidth',2);
p1=plot(34.57,3.981,'bo','MarkerSize', 6);
p2=plot(70.64,3.984,'ro','MarkerSize', 6);
p3=plot(82.67,3.98,'go','MarkerSize', 6);
legend([p1 p2 p3],{'34.6 ka' '70.6 ka' '82.7 ka'},'Location','E')
xlim([0 250])
xlabel('time [ka]')
ylabel('thickness [km]')
title('Equatorial Ice Thickness')

subplot(1,2,2)
hold on 
grid on

%p11=plot(t12./10^3,Hd12(1:100:end,1)./10^3,'b-','LineWidth',2);
p11=plot(t13./10^3,H12b(:,1)./10^3,'b-','LineWidth',2);
p22=plot(t13./10^3,Hd13(1:500:end,1)./10^3,'r-','LineWidth',2);
p33=plot(t13./10^3,H14i(:,1)./10^3,'g-','LineWidth',2);
legend([p11 p22 p33], {'\eta_0=10^{12}','\eta_0=10^{13}','\eta_0=10^{14}'},'Location','NE')
xlabel('time [ka]')
xlim([0 250])
ylabel('thickness [km]')
title('Polar Ice Thickness')
