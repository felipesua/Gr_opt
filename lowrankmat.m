cd 'C:\Users\Felipe\Downloads\projects\Gr'

%% 1 - Rank 1 
clc, clear, close all

x = rand(3,1); x = x/norm(x);
x = [0;0;1];
y = rand (3,1); y = y/norm(y);

[X,Y,Z] = sphere(32);
T = 1-(X*x(1)+Y*x(2)+Z*x(3)).^2;


%% Plotting
f1= figure(1); subplot(111), clc
surf(X,Y,Z,T)

xlabel('$\textbf{x}_1$'), ylabel('$\textbf{x}_2$'), zlabel('$\textbf{x}_3$'), %title('$f(\textbf{x})$')
pbaspect([1 1 1]); axis equal, shading interp, colormap parula;
cb = colorbar; set(cb,'position',[.89 .3 ,.04 .35]) 
set(f1,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9+2^5, 2^9])
set(gca,'fontsize', 20)

%% Persistence
export_fig './figures/rank1.pdf' -pdf -transparent -c[0,Inf,Inf,Inf]
saveas(f1,'./figures/rank1.fig')

%% 2 - Full obs landscape
clc, clear

m = 100;
n = 200;
r = 6;
tiks = 32;
p = .9; 

U = orth(randn(m,r));
V = orth(randn(m,r));  [U,V] = orthogonal_project_normalize(U,V);
M = U*diag([1 1 1 1.4 1.4 1.4])*orth(rand(n,r))';

O = rand(m,n) < p;
[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));
f = partial_projectionDist(tiks,U,V,M,O,p);

%% Plotting
f2 = figure(1); subplot(111), clc
surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, hold on
xlabel('$\theta_1$'), ylabel('$\theta_2$')

pbaspect([1 1 1]), cb = colorbar; set(cb,'position',[.875 .3 ,.04 .425]), caxis([0,4.5])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'})
set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), 
set(f2,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9+2^7, 2^9])
set(gca,'fontsize', 20)
hold off


%% Persistence
export_fig './figures/full_landscape.pdf' -pdf -transparent -c[0,Inf,Inf,Inf]
saveas(f2,'./figures/full_landscape.fig')

%% 2 Landscape/geodesics
clc, clear

m = 100;
n = 200;
r = 6;
tiks = 32;
pts = 20;
p = 1; 

U = orth(randn(m,r));
V = orth(randn(m,r));  [U,V] = orthogonal_project_normalize(U,V);
M = U*diag([1 .5 1 .5 1 .5])*orth(rand(n,r))';

O = rand(m,n) < p;

f = partial_projectionDist(tiks,U,V,M,O,p);
[t1,t2,f1] = geodesic2([4/8*pi/2,0],[0,5/8*pi/2],V,U,O,M,pts,p);
[t3,t4,f2] = geodesic2([0,2/8*pi/2],[2/4*pi/2,3/4*pi/2],V,U,O,M,pts,p);
[t5,t6,f3] = geodesic2([3/8*pi/2,0],[pi/2,pi/2],V,U,O,M,pts,p);
[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));

%% Plotting
fig3 = figure(1); clc, 
fonts = 22;

subplot(2,3,[1 2 4 5])
surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, hold on
plot3(t1,t2,f1+.01,'.-','Markersize',8,'color', 'magenta')
plot3(t3,t4,f2+.01,'.-','Markersize',8,'color', [.54 0 .54])
plot3(t5,t6,f3+.01,'.-','Markersize',8,'color', [.85 .125 .125])
plot3(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',24,'color', 'magenta')
plot3(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),f2(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',24,'color', [.54 0 .54])
plot3(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),f3(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',24,'color', [.85 .125 .125])

pbaspect([1 1 .85]), cb = colorbar; set(cb,'position',[.075 .3 ,.025 .45]), caxis([0,2])
set(gca,'ztick',0:1:4,'zticklabel',{'0','1','2','3','4'})
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'})
set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
view(-40,22), xlabel('$\theta_1$'), ylabel('$\theta_2$')
set(gca,'fontsize', fonts), pbaspect([1 1 1]), 
hold off


subplot(2,3,3)
imagesc(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks),f), hold on,  caxis([0,2])
plot(t1,t2,'.-','Markersize',8,'color', 'magenta')
plot(t3,t4,'.-','Markersize',8,'color', [.54 0 .54])
plot(t5,t6,'.-','Markersize',8,'color', [.85 .125 .125])
plot(t1(find(t2<pi/4 & t1<pi/4, 1,'last')),t2(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-','Markersize',24,'color', 'magenta')
plot(t3(find(t3<pi/4 & t4<pi/4, 1,'last')),t4(find(t3<pi/4 & t4<pi/4, 1,'last')),'.-','Markersize',24,'color', [.54 0 .54])
plot(t5(find(t5<pi/4 & t6<pi/4, 1,'last')),t6(find(t5<pi/4 & t6<pi/4, 1,'last')),'.-','Markersize',24,'color', [.85 .125 .125])

pbaspect([1 1 1])
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'})
set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
set(gca,'YDir','normal'), set(gca,'xaxisLocation','bottom')
xlabel('$\theta_1$'), ylabel('$\theta_2$')
set(gca,'fontsize', fonts)
hold off


subplot(2,3,6)
plot(linspace(0,1,pts),flip(f1),'.-','color','magenta'), hold on
plot(linspace(0,1,pts),flip(f2),'.-','color',[.54 0 .54]), 
plot(linspace(0,1,pts),flip(f3),'.-','color',[.85 .125 .125]), 
plot(1-(find(t2<pi/4 & t1<pi/4, 1,'last')-1)/(pts-1),f1(find(t2<pi/4 & t1<pi/4, 1,'last')),'.-', 'Markersize',20, 'color','magenta'), 
plot(1-(find(t4<pi/4 & t3<pi/4, 1,'last')-1)/(pts-1),f2(find(t4<pi/4 & t3<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.54 0 .54]), 
plot(1-(find(t6<pi/4 & t5<pi/4, 1,'last')-1)/(pts-1),f3(find(t6<pi/4 & t5<pi/4, 1,'last')),'.-', 'Markersize',20, 'color',[.85 .125 .125]), 

pbaspect([1 1 1])
grid on,, xlabel('$t$'),ylabel('$f_{\textbf{full}}(X_t)$')
set(gca,'fontsize', fonts)
hold off

set(fig3,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9+2^9, 2^9])

%% Persistence
export_fig './figures/full_geodesics.pdf' -pdf -transparent -c[Inf,Inf,1,Inf]
saveas(fig3,'./figures/full_geodesics.fig')

%% 3 Rate
clc, clear

m = 100;
n = 200;
r = 2;
tiks = 32;
pts = 30;
eta = .2;

Um = incohere(m,r,5000); U =Um;
% U = orth(randn(m,r));
V = orth(randn(m,r));  [U,V] = orthogonal_project_normalize(U,V);

M = U*diag([1 1.2])*orth(rand(n,r))';
f = projectionDist(tiks,U,V,M);

[dX,dY] = meshgrid(linspace(pi/2,0,tiks),linspace(pi/2,0,tiks));

%% Plotting
fig4 = figure(1); clc

subplot(3,2,1:4)
surf(dX,dY,f,'EdgeColor',[.75 .75 .75]), colormap parula, hold on, 
pbaspect([1 1 1]), view(-40,22)
set(gca,'xtick',0:pi/4:pi/2,'xticklabel',{'0','\pi/4','\pi/2'})
set(gca,'ytick',0:pi/4:pi/2,'yticklabel',{'0','\pi/4','\pi/2'})
xlabel('$\theta_1$'), ylabel('$\theta_2$')
cb = colorbar; set(cb,'position',[.89 .52 ,.035 .26]), caxis([0,1.5])
set(gca,'fontsize', 14)


color = 0;
for i = linspace(0.1,pi/2-.1,3)
s0 = .95*cos(i)^.35*pi/2; t0 = .95*sin(i)^.15*pi/2;
X = U*diag([cos(s0) cos(t0)]) + V*diag([sin(s0) sin(t0)]);
[ss,tt,ff,fp,inc] = proj_trajectory(X,U,V,M,eta,pts);

subplot(3,2,1:4)
if color == 0
    color = 'magenta';
elseif strcmp(color,'magenta')
    color = [.54 0 .54];
else
    color = [.85 .125 .125];
end
plot3(ss,tt,fp+.01,'.-','Markersize',14,'color',color)

subplot(3,2,5)
semilogy(ff, '.-','color', color, 'Markersize',5),  hold on, grid on

subplot(3,2,6)
plot(10*sqrt(inc), '.-','color',color), hold on

end

subplot(3,2,1:4), hold off
text(-.5,-.5,0, '(a)','fontsize',15)

subplot(3,2,5), 
text(11,3e-6, '(b)','fontsize',15), xlabel('$k$'), ylabel('$\log f(X_k)$')
set(gca,'fontsize', 13), ylim([10^(-4),1])
grid on, hold off

subplot(3,2,6)
text(11,1.45, '(c)','fontsize',15), xlabel('$k$'), ylabel('$10\cdot\max_i ||e_i^\top X_k||$')
set(gca,'fontsize', 13)
grid on, hold off

set(fig4,'WindowStyle','normal')
set(gcf, 'Position', [2^3, 2^8, 2^9, 2^9+1.75*2^7])

%%
export_fig './figures/GD-rate-incoherence.pdf' -pdf -transparent -c[Inf,Inf,10,Inf]
saveas(fig4,'./figures/GD-rate-incoherence.fig')

%% Aux

function [X,Y,Z] = pos_sphere(n)
[V,U] = meshgrid(linspace(0,pi/2,n), linspace(0,pi/2,n));

X = cos(V) .* sin(U);
Y =  sin(V);
Z = cos(V) .* cos(U);
end


% [X,Y,Z] = pos_torus(.1,1,tiks); surf(X,Y,Z,f)
% xlabel('$\theta_1$'), ylabel('$\theta_2$'), title('$f(\textbf{x})$')
% pbaspect([1 1 1]), axis equal,  %shading interp, colorbar 
function [X,Y,Z] = pos_torus(R,r,n)

[V,U] = meshgrid(-pi/4+linspace(0,pi/2,n), 0*pi/4+linspace(0,pi/2,n));

X = (R+r*cos(V)).* cos(U);
Y = (R+r*cos(V)).*sin(U);
Z = r*sin(V);

end

function f = partial_Boumal(n,U,V,M,O,p,l)
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
%         f(si,ti) = 1/p*.5*norm(O.*(M - X*X'*M),'fro')^2; % TODO correct this
        A = ((1-l^2)* kron(eye(m), X')*diag(O(:))*kron(eye(m),X) + l^2*eye(m*r));
        Y = (M.*O)'*U;
        Y = reshape(A\reshape(Y,[],1),m,r);
        f(si,ti) = .5/p*norm(O.*(X*Y' - O.*M),'fro')^2 + .5/p*l^2*(norm(Y,'fro')^2-norm(O.*(X*Y'),'fro')^2);
    end
end
end

function f = partial_projectionDist(n,U,V,M,O,p)
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
%         X(:,1) = V(:,1)*cos(t) + U(:,1)*sin(t);
%         X(:,2) = V(:,2)*cos(s) + U(:,2)*sin(s);
        X = V*diag([cos(t)*ones(1,r) cos(s)*ones(1,r)]) + ...
            U*diag([sin(t)*ones(1,r) sin(s)*ones(1,r)]);
        f(si,ti) = 1/p*.5*norm(O.*(M - X*X'*M),'fro')^2;
    end
end
end


function f = projectionDist(n,U,V,M)
f = zeros(n,n);

X = 0*U;
ss = linspace(0,pi/2,n);
tt = linspace(0,pi/2,n);
for si = 1:n
    for ti = 1:n
        s = ss(si);
        t = tt(ti);
        X(:,1) = V(:,1)*cos(t) + U(:,1)*sin(t);
        X(:,2) = V(:,2)*cos(s) + U(:,2)*sin(s);
        f(si,ti) = .5*norm(M - X*X'*M,'fro')^2;
    end
end

end

function [t1,t2,f] = geodesic2(a,b,V,U,O,M,tiks,p)
t1 = linspace(a(1),b(1),tiks);
t2 = linspace(a(2),b(2),tiks);
r = size(V,2)/2;
f = 0*t1;
for i = 1:tiks
    X = U*diag([cos(t1(i))*ones(1,r) cos(t2(i))*ones(1,r)]) +...
        V*diag([sin(t1(i))*ones(1,r) sin(t2(i))*ones(1,r)]);
    f(i) = 1/p*.5*norm(O.*(M-X*X'*M),'fro')^2;
end
end

function [t1,t2,f] = geodesic(s1,s2,V,U,O,M,tiks,p)
t1 = linspace(0,s1,tiks);
t2 = linspace(0,s2,tiks);
r = size(V,2)/2;
f = 0*t1;
for i = 1:tiks
    X = U*diag([cos(t1(i))*ones(1,r) cos(t2(i))*ones(1,r)]) +...
        V*diag([sin(t1(i))*ones(1,r) sin(t2(i))*ones(1,r)]);
    f(i) = 1/p*.5*norm(O.*(M-X*X'*M),'fro')^2;
end
end


function [s,t] = rank2Project(X,V,U)
s = zeros(1,2);
for i = 1:2
u = X'*V(:,i);
v = X'*U(:,i);

a = norm(u)^2;
b = 2*u'*v;
c = norm(v)^2;

if b<=0
    if a>=c
        s(i) = 0;
    else
        s(i) = pi/2;
    end
else
    s(i) = .5*(pi/2 - atan((a-c)/b));
end
end
t = s(2);
s = s(1);
end

function [ss,tt,f,fp,inc] = proj_trajectory(X0,U,V,M,eta,iter)
% TODO add incoherence marker
ss = zeros(iter,1);
tt = zeros(iter,1);
f = zeros(iter,1);
fp = zeros(iter,1);
X = X0;
m = size(M,1);
inc = zeros(iter,1);
for i = 1:iter
%  X
grad = (eye(m) - X*X')*M*M'*X;
X = orth(X + eta*grad);
f(i) = .5*norm(M-X*X'*M,'fro')^2;
inc(i) = max(sum(X.^2,2));

% pX
[s,t] = rank2Project(X,U,V);
ss(i) = s; tt(i) = t;
pX = U*diag([cos(s) cos(t)]) + V*diag([sin(s) sin(t)]);
fp(i) = .5*norm(M-pX*pX'*M,'fro')^2;

end
end

function [U,X] = orthogonal_project_normalize(U,X)
[A,S,B] = svd(X' * U);
X = X*A;
U = U*B;
X = X - U*U'*X;
X = X./sqrt(sum(X.^2,1));
end

function Um = incohere(m,r,iter)
inc = 10;
Um = orth(rand(m,r));
for i = 1:iter
    U = orth(randn(m,r));
    if  max(sum(U.^2,2))*m/r < inc
        Um = U;
        inc = max(sum(U.^2,2))*m/r;
    end
end
end

%% Gradient descent 
% m = 25;
% r = 4;
% r0 = 0;
% eta = 1e-3;
% iter = 1e4;
% res = ones(1,iter);
% res2 = ones(1,iter);
% % sig = .1;
% % rng(823)
% 
% % U = randn(m,r);
% U = [eye(r); zeros(m-r,r)];
% 
% % M = U*diag(1+rand(1,r))*orth(rand(m,r))';
% % M = U*diag(ones(1,r))*orth(rand(m,r))';
% M = U*diag(linspace(1,2,r))*orth(rand(m,r))';
% 
% X = orth(randn(m,r));
% X(1:m,1:r0) = [eye(r0); zeros(m-r0,r0)];
% 
% %
% for i = 1:iter
% grad = (eye(m) - X*X')*M*M'*X;
% X = orth(X + eta*grad);
% res(i) = norm(X*X'*M-M,'fro');
% res2(i) =  norm(grad,2);
% end
% 
% figure(1)
% yyaxis left,  plot(res), grid on,   ylim([0,3]), xlim([1,iter])
% yyaxis right,  plot(res2), grid on
% 

