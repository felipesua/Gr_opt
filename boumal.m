cd 'C:\Users\Felipe\Downloads\projects\Gr'
addpath '.'/shaded_plots/
%% 0 Partial concetration tests: col by col concentration
% clc, clear

m = 5e2;
n = 5e2;
r = 8;
tiks = 20;

U = orth(randn(m,r));
U = incoherentU(m,r);

V = orth(randn(m,r));  [U,V] = orthogonal_project_normalize(U,V);
M = U*diag(linspace(1,2,r))*orth(rand(n,r))';

ffull = partial_col(tiks,U,V,M,true(m,n),1,1:n)+1e-4;
%% Random mask
p = .5; %0.026 for scenario 2, plot each scenario against full obs with geodesics;
O = rand(m,n) < p;

[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));
f = partial_col(tiks,U,V,M.*O,O,p,1:n)+1e-4;

% %%
subplot(221),  surf(dX,dY,ffull, 'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22)
subplot(222),  surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22)
subplot(223), surf(dX,dY,ffull./f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22), zlim([.75,1.25])

%% Principal angle concentrationo
clc
% U = orth(randn(m,r));
% X = orth(randn(m,r));
U = incoherentU(m,r);
X = incoherentU(m,r);
[U,X] = orthogonal_project_normalize(U,X);
ths = [ones(1,r/2)*1*pi/8 linspace(0*pi/8,4*pi/8,r/2)];
% ths = ones(1,r) * pi/18;
X = U*diag(cos(ths)) + X*diag(sin(ths));

O = rand(m,n) < .15;
reps = 5e2;

svd_t = svd(X' * U);
sn = sum(svd_t);

% %%
svd_c = zeros(r,reps);
sn_ = zeros(1,reps);
for i = 1:reps
Y = X(O(:,i),:);
svd_c(:,i) = svd(pinv(Y' * Y) * Y' * U(O(:,i),:));
sn_(i) = sum(svd_c(:,i));
end

[svd_t, mean(svd_c,2), abs(svd_t - mean(svd_c,2)), std(svd_c')']

[sn, mean(sn_) std(sn_)]
subplot(224), plot(abs(svd_t - mean(svd_c,2))./svd_t, '-o'), grid on,  ylim([0,.050]), xlim([1 r])










%% (Everything below is OLD) 6 - Scenario Concentration
clc, clear
m = 1e3;
n = 1e3;
r = 10;

% M and U 
U = orth(randn(m,r)); 
% M = U* diag( linspace(1,2,r) ) *orth(rand(n,r))'; % Well-conditioned
M = U* diag( sqrt(m*n)*exp(linspace(0,-5,r)) ) *orth(rand(n,r))'; % Bad-conditioned

% Random pts in Gr
X1 = orth(randn(m,r)); [U,X1] = orthogonal_project_normalize(U,X1); X1 = rotate_angles(U,X1, pi/2*[ones(1,r/2) .5*ones(1,r/2)] );
X2 = orth(randn(m,r)); [U,X2] = orthogonal_project_normalize(U,X2); X2 = rotate_angles(U,X2, pi/2*ones(1,r));
X3 = orth(randn(m,r)); [U,X3] = orthogonal_project_normalize(U,X3); X3 = rotate_angles(U,X3, pi/4*ones(1,r));
X4 = orth(randn(m,r)); [U,X4] = orthogonal_project_normalize(U,X4); X4 = rotate_angles(U,X4, pi/2*[3/4*ones(1,1) 1/4*ones(1,r-1)] );
X5 = orth(randn(m,r)); [U,X5] = orthogonal_project_normalize(U,X5); X5 = rotate_angles(U,X5, pi/2*[3/4*ones(1,3) 1/4*ones(1,r-3)] );
X6 = orth(randn(m,r)); [U,X6] = orthogonal_project_normalize(U,X6); X6 = rotate_angles(U,X6, pi/2*[3/4*ones(1,5) 1/4*ones(1,r-5)] );
X7 = orth(randn(m,r)); [U,X7] = orthogonal_project_normalize(U,X7); X7 = rotate_angles(U,X7, pi/2*[3/4*ones(1,7) 1/4*ones(1,r-7)] );
X8 = orth(randn(m,r)); [U,X8] = orthogonal_project_normalize(U,X8); X8 = rotate_angles(U,X8, pi/2*[3/4*ones(1,9) 1/4*ones(1,r-9)] );
X9 = orth(randn(m,r)); 

X = {X1, X2, X3, X4, X5, X6, X7, X8, X9};
clearvars X1 X2 X3 X4 X5 X6 X7 X8 X9 

% Trial parameters
ps = 10.^[linspace(-3,-1,9) log(.5)/log(10)];

% results
trials = 1e2;
fs = zeros(10,trials,9);


for i = 1:10
    for j = 1:trials
        O = rand(m,n) < ps(i);
        for k = 1:9
            if i == 1 && j == 1 && k ==1
                tic
                Y = partial_inverse(X{k},O,M.*O);
                fs(i,j,k) = .5/ps(i)*norm(O.*(X{k}*Y' - M),'fro')^2;
                toc * 10* 9 * trials/60
            else
                Y = partial_inverse(X{k},O,M.*O);
                fs(i,j,k) = .5/ps(i)*norm(O.*(X{k}*Y' - M),'fro')^2;
            end
        end
    end
end

%% Plotting
fig6 = figure(1); clc
fs = fs/(m*n);

p1 = semilogx(ps, mean(fs(:,:,4)')','o-','LineWidth',1,'MarkerSize',2); hold on;
p2 = semilogx(ps, mean(fs(:,:,4)')'+1.96*std(fs(:,:,4)')','-','LineWidth',2); p2.Color = p1.Color; p2.Color(4) = .5;
p2 = semilogx(ps, mean(fs(:,:,4)')'-1.96*std(fs(:,:,4)')','-','LineWidth',2); p2.Color = p1.Color; p2.Color(4) = .5;
text(5e-1, mean(fs(end,:,4)')','$X_1$','horizontal','left', 'vertical','middle', 'FontSize', 13)

for i = 5:9
p1 = semilogx(ps, mean(fs(:,:,i)')','o-','LineWidth',1,'MarkerSize',2); 
p2 = semilogx(ps, mean(fs(:,:,i)')'+1.96*std(fs(:,:,i)')','-','LineWidth',2); p2.Color = p1.Color; p2.Color(4) = .5;
p2 = semilogx(ps, mean(fs(:,:,i)')'-1.96*std(fs(:,:,i)')','-','LineWidth',2); p2.Color = p1.Color; p2.Color(4) = .5;
text(5e-1, mean(fs(end,:,i)')'+.4e2*(i==7)-.3e2*(i==5),strcat('$X_',num2str(i-3),'$'),'horizontal','left', 'vertical','middle', 'FontSize', 13)
end
grid on, xlabel('$p$'), ylabel('$f(\textbf{X})$'),%ylim([0,8])
set(gca,'fontsize', 15),% pbaspect([1 1 1]), 
hold off

set(fig6,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9, 2^8+2^6])

%% Persistence
cd 'C:\Users\Felipe\Downloads\projects\Gr'
export_fig './figures/sc3_f_concentration.pdf' -pdf -transparent -c[0,Inf,Inf,Inf]
saveas(fig6,'./figures/sc3_f_concentration.fig')

%% 5 - Partial Landscape
clc, clear

m = 1e2;
n = 1e2;
r = 4;
tiks = 20;

U = orth(randn(m,r));
V = orth(randn(m,r));  [U,V] = orthogonal_project_normalize(U,V);
M = U*diag(linspace(1,2,r))*orth(rand(n,r))';

ffull = partial_Boumal0(tiks,U,V,M,ones(m,n),1);
%%
p = .5; %0.026 for scenario 2, plot each scenario against full obs with geodesics;
O = rand(m,n) < p;

[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));
f = partial_Boumal0(tiks,U,V,M.*O,O,p);

%%
subplot(221),  surf(dX,dY,ffull, 'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22)
subplot(222),  surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22)
%%

steps = 16;
[t1,t2,f1] = geodesic([4/8*pi/2,0],[0,5/8*pi/2],V,U,O,O.*M,steps,p);
[t3,t4,f2] = geodesic([0,2/8*pi/2],[2/4*pi/2,3/4*pi/2],V,U,O,O.*M,steps,p);
[t5,t6,f3] = geodesic([3/8*pi/2,0],[pi/2,pi/2],V,U,O,O.*M,steps,p);


%%
p1.f = f;
p1.f1 = f1;
p1.f2 = f2;
p1.f3 = f3;

%% Plotting
fig5 = figure(1); clc

tiks = 16; steps = 16;
[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));
t1 = linspace(4/8*pi/2,0,tiks); t2 = linspace(0,5/8*pi/2,tiks);
t3 = linspace(0*pi/2,2/4*pi/2,tiks); t4 = linspace(2/8*pi/2,3/4*pi/2,tiks);
t5 = linspace(3/8*pi/2,pi/2,tiks); t6 = linspace(0,pi/2,tiks);

p = p10;
q = p1;
offs = 1e-2;

p.f = p.f*1e-6; p.f1 = p.f1*1e-6; p.f2 = p.f2*1e-6; p.f3 = p.f3*1e-6; % only for scenario 3
q.f = q.f*1e-6; q.f1 = q.f1*1e-6; q.f2 = q.f2*1e-6; q.f3 = q.f3*1e-6;


subplot(221),  surf(dX,dY,q.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,q.f1+offs,':','linewidth',1.5,'color', 'magenta')
plot3(t3,t4,q.f2+offs,':','linewidth',1.5,'color', [.54 0 .54])
plot3(t5,t6,q.f3+offs,':','linewidth',1.5,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),q.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),q.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),q.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$'), 
set(gca,'fontsize', 13), pbaspect([1 1 1]),% zlim([0,12])
hold off
%%

subplot(222),  surf(dX,dY,p.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,p.f1+offs,'.-','Markersize',8,'color', 'magenta')
plot3(t3,t4,p.f2+offs,'.-','Markersize',8,'color', [.54 0 .54])
plot3(t5,t6,p.f3+offs,'.-','Markersize',8,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),p.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),p.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),p.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$')
set(gca,'fontsize', 13), pbaspect([1 1 1]), %zlim([0,12])
hold off

% subplot(2,2,[3,4])
subplot(2,2,3)
plot(linspace(0,1,steps),flip(q.f1),':','linewidth',2,'color','magenta'), hold on
plot(linspace(0,1,steps),flip(q.f2),':','linewidth',2,'color',[.54 0 .54]), 
plot(linspace(0,1,steps),flip(q.f3),':','linewidth',2,'color',[.85 .125 .125]), 
plot(1-(find(t2<pi/4 & t1<pi/4, 1,'last')-1)/(steps-1),q.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-', 'Markersize',20, 'color','magenta'), 
plot(1-(find(t4<pi/4 & t3<pi/4, 1,'last')-1)/(steps-1),q.f2(find(t4<pi/4 & t3<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.54 0 .54]), 
plot(1-(find(t6<pi/4 & t5<pi/4, 1,'last')-1)/(steps-1),q.f3(find(t6<pi/4 & t5<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.85 .125 .125]), 
grid on,, xlabel('$t$'),ylabel('$f(\textbf{X}_t)$'), %ylim([0,12])
set(gca,'fontsize', 13), %pbaspect([1 1 1]), 
hold off

subplot(2,2,4)
plot(linspace(0,1,steps),flip(p.f1),'.-','color','magenta'), hold on
plot(linspace(0,1,steps),flip(p.f2),'.-','color',[.54 0 .54]), 
plot(linspace(0,1,steps),flip(p.f3),'.-','color',[.85 .125 .125]), 
plot(1-(find(t2<pi/4 & t1<pi/4, 1,'last')-1)/(steps-1),p.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-', 'Markersize',20, 'color','magenta'), 
plot(1-(find(t4<pi/4 & t3<pi/4, 1,'last')-1)/(steps-1),p.f2(find(t4<pi/4 & t3<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.54 0 .54]), 
plot(1-(find(t6<pi/4 & t5<pi/4, 1,'last')-1)/(steps-1),p.f3(find(t6<pi/4 & t5<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.85 .125 .125]), 
grid on,, xlabel('$t$'),ylabel('$f_{\textbf{full}}(\textbf{X}_t)$'), %ylim([0,12])
set(gca,'fontsize', 13), %pbaspect([1 1 1]), 
hold off

set(fig5,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9+2^7, 2^9])

%% Persistence 

export_fig './figures/sc3-full-boumal.pdf' -pdf -transparent -c[0,Inf,Inf,Inf]
% export_fig './figures/sc1-full-boumal.eps' -eps -transparent -c[0,Inf,Inf,Inf]
% saveas(fig5,'./figures/sc1-full-boumal.fig')


%% 7 - Landscape varying p

fig7 = figure(1); clc, 

tiks = 16; steps = 16;
[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));
t1 = linspace(4/8*pi/2,0,tiks); t2 = linspace(0,5/8*pi/2,tiks);
t3 = linspace(0*pi/2,2/4*pi/2,tiks); t4 = linspace(2/8*pi/2,3/4*pi/2,tiks);
t5 = linspace(3/8*pi/2,pi/2,tiks); t6 = linspace(0,pi/2,tiks);

p = p10;
q = p1;
r = p01;
s = p001;
offs = 1e-2;

% p.f = p.f*1e-6; p.f1 = p.f1*1e-6; p.f2 = p.f2*1e-6; p.f3 = p.f3*1e-6;
% q.f = q.f*1e-6; q.f1 = q.f1*1e-6; q.f2 = q.f2*1e-6; q.f3 = q.f3*1e-6;
% r.f = r.f*1e-6; r.f1 = r.f1*1e-6; r.f2 = r.f2*1e-6; r.f3 = r.f3*1e-6;
% s.f = s.f*1e-6; s.f1 = s.f1*1e-6; s.f2 = s.f2*1e-6; s.f3 = s.f3*1e-6;

subplot(221),  surf(dX,dY,p.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,p.f1+offs,':','linewidth',1.5,'color', 'magenta')
plot3(t3,t4,p.f2+offs,':','linewidth',1.5,'color', [.54 0 .54])
plot3(t5,t6,p.f3+offs,':','linewidth',1.5,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),p.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),p.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),p.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$'), 
set(gca,'fontsize', 13), pbaspect([1 1 1]),% zlim([0,12])
hold off


subplot(222),  surf(dX,dY,q.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,q.f1+offs,':','linewidth',1.5,'color', 'magenta')
plot3(t3,t4,q.f2+offs,':','linewidth',1.5,'color', [.54 0 .54])
plot3(t5,t6,q.f3+offs,':','linewidth',1.5,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),q.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),q.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),q.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$'), 
set(gca,'fontsize', 13), pbaspect([1 1 1]),% zlim([0,12])
hold off

subplot(223),  surf(dX,dY,r.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,r.f1+offs,':','linewidth',1.5,'color', 'magenta')
plot3(t3,t4,r.f2+offs,':','linewidth',1.5,'color', [.54 0 .54])
plot3(t5,t6,r.f3+offs,':','linewidth',1.5,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),r.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),r.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),r.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$'), 
set(gca,'fontsize', 13), pbaspect([1 1 1]),% zlim([0,12])
hold off

subplot(224),  surf(dX,dY,s.f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),   hold on
plot3(t1,t2,s.f1+offs,':','linewidth',1.5,'color', 'magenta')
plot3(t3,t4,s.f2+offs,':','linewidth',1.5,'color', [.54 0 .54])
plot3(t5,t6,s.f3+offs,':','linewidth',1.5,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),s.f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',20,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),s.f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',20,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),s.f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',20,'color', [.85 .125 .125])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$'), 
set(gca,'fontsize', 13), pbaspect([1 1 1]),% zlim([0,12])
hold off

set(fig7,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9+2^6, 2^9])
%%

export_fig './figures/sc1-varp.pdf' -pdf -transparent -c[0,Inf,Inf,Inf]
saveas(fig7,'./figures/sc3-varp.fig')


%% Auxiliary functions

function f = partial_col(n,U,V,M,O,p,j)
f = zeros(n,n);
m = size(U,1);
r = size(V,2)/2;

X = 0*U;
ss = linspace(0,pi/2,n);
tt = linspace(0,pi/2,n);
for si = 1:n
    for ti = 1:n
        s = ss(si);
        t = tt(ti);
        X = V*diag([cos(t)*ones(1,r) cos(s)*ones(1,r)]) + ...
            U*diag([sin(t)*ones(1,r) sin(s)*ones(1,r)]);
%         A = kron(eye(m), X')*diag(O(:))*kron(eye(m),X);
%         Y = (M.*O)'*X;
%         Y = reshape(pinv(A)*reshape(Y',[],1),2*r,m)';
        Y = partial_inverse(X,O,M);
        res = O.*(X*Y' - M);
        res = res(:,j);
        f(si,ti) = .5/p*norm(res,'fro')^2;
    end
end
end

function [t1,t2,f] = geodesic(a,b,V,U,O,M,steps,p)
t1 = linspace(a(1),b(1),steps);
t2 = linspace(a(2),b(2),steps);
r = size(V,2)/2;
f = 0*t1;
for i = 1:steps
    X = U*diag([cos(t1(i))*ones(1,r) cos(t2(i))*ones(1,r)]) +...
        V*diag([sin(t1(i))*ones(1,r) sin(t2(i))*ones(1,r)]);
    Y = partial_inverse(X,O,M);
    f(i) = .5/p*norm(O.*(X*Y' - M),'fro')^2;
end
end

function f = partial_Boumal0(n,U,V,M,O,p)
f = zeros(n,n);
m = size(U,1);
r = size(V,2)/2;

X = 0*U;
ss = linspace(0,pi/2,n);
tt = linspace(0,pi/2,n);
for si = 1:n
    for ti = 1:n
        s = ss(si);
        t = tt(ti);
        X = V*diag([cos(t)*ones(1,r) cos(s)*ones(1,r)]) + ...
            U*diag([sin(t)*ones(1,r) sin(s)*ones(1,r)]);
%         A = kron(eye(m), X')*diag(O(:))*kron(eye(m),X);
%         Y = (M.*O)'*X;
%         Y = reshape(pinv(A)*reshape(Y',[],1),2*r,m)';
        Y = partial_inverse(X,O,M);
        res = O.*(X*Y' - M);
        res = res(:,1);
        f(si,ti) = .5/p*norm(res,'fro')^2;
    end
end
end

function Y = partial_inverse(X,O,M)
    m = size(O,2);
    Y = (M)'*X;
    for i = 1:m
        Xp = X(O(:,i),:);
        Y(i,:) = Y(i,:)*pinv(Xp'*Xp);
%         Y(i,:) = Y(i,:)*pinv(X'*diag(O(:,i))*X);
    end
end


function f = partial_Boumal_old(n,U,V,M,O,p,l)
f = zeros(n,n);
m = size(U,1);
r = size(V,2)/2;

X = 0*U;
ss = linspace(0,pi/2,n);
tt = linspace(0,pi/2,n);
for si = 1:n
    for ti = 1:n
        s = ss(si);
        t = tt(ti);
        X = V*diag([cos(t)*ones(1,r) cos(s)*ones(1,r)]) + ...
            U*diag([sin(t)*ones(1,r) sin(s)*ones(1,r)]);
        A = ((1-l^2)* kron(eye(m), X')*diag(O(:))*kron(eye(m),X) + l^2*eye(m*2*r));
        Y = (M.*O)'*X;
        if l == 0
            Y = reshape(pinv(A)*reshape(Y',[],1),2*r,m)';
        else
            Y = reshape(A\reshape(Y',[],1),2*r,m)';
        end
        f(si,ti) = .5/p*norm(O.*(X*Y' - O.*M),'fro')^2 + .5/p*l^2*(norm(Y,'fro')^2-norm(O.*(X*Y'),'fro')^2);
    end
end
end

% f = proxy(tiks,U,V,M,O,p);
% subplot(121),  surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, view(-40,22),  
% title(strcat('Proxy: $p=', num2str(p), '$' ))
% zlim([0,r]), set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'}), set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})

function f = proxy(n,U,V,M,O,p)
f = zeros(n,n);
m = size(U,1);
r = size(V,2)/2;

X = 0*U;
ss = linspace(0,pi/2,n);
tt = linspace(0,pi/2,n);
for si = 1:n
    for ti = 1:n
        s = ss(si);
        t = tt(ti);
        X = V*diag([cos(t)*ones(1,r) cos(s)*ones(1,r)]) + ...
            U*diag([sin(t)*ones(1,r) sin(s)*ones(1,r)]);
        Y = M'*X;
        f(si,ti) = 1/p*.5*norm(O.*(X*Y' - M),'fro')^2;
    end
end
end


function Um = incoherentU(m,r)
inc = Inf;
Um = orth(rand(m,r));
for i = 1:500
    U = orth(randn(m,r));
    if  max(sum(U.^2,2))*m/r < inc
        Um = U;
        inc = max(sum(U.^2,2))*m/r;
    end
end
end

function Um = coherentU(m,r)
inc = -Inf;
Um = orth(rand(m,r));
for i = 1:500
    U = orth(randn(m,r));
    if  min(sum(U.^2,2))*m/r > inc
        Um = U;
        inc = min(sum(U.^2,2))*m/r;
    end
end
end

function [U,X] = orthogonal_project_normalize(U,X)
[A,S,B] = svd(X' * U);
X = X*A;
U = U*B;
X = X - U*U'*X;
X = X./sqrt(sum(X.^2,1));
end

function X = rotate_angles(U,X,angles)
X = U*diag(cos(angles)) + X*diag(sin(angles));
end