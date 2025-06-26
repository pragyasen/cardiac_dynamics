%figure1
%graphing v-nullcline: dv/dt = 0
fimplicit(@(va,ha) (-8*va)*(va-0.15)*(va-1)-(va*ha),'b')
hold on
%graphing h-nullclines: dh/dt = 0
fimplicit(@(vb,hb) (0.002 +(0.2 * hb/(vb + 0.3))),'r')
hold on
fimplicit(@(vc,hc) (-8*vc.^2 + 9.2*vc - hc), 'r')
hold on
%setting axis boundaries:
axis([-0.05 1.2 -0.05 2.75])
xlabel('v')
ylabel('h')
legend('v-nullcline','h-nullcline');
%plotting equilibrium point
plot(0,0,'.black','MarkerSize', 20)
text(0.01,0.1,'(0,0)');

%figure2
close all;
clear all;
clc;
vmin = -0.05; vmax = 1.2;
hmin = -0.05; hmax = 2.75;
%step size for v and h
vstep = 0.5; hstep = 0.5;
%generating mesh for plotting
[v, h] = meshgrid(vmin:0.1:vmax, hmin:0.1:hmax);
%setting parameter values
a = 0.15; k = 8; ep = 0.002; m1 = 0.2; m2 = 0.3;
%system of equations dv/dt = f(v,h) and dh/dt = g(v,h)
dv = -k*v.*(v-a).*(v-1) - v.*h;
dh = (ep + (m1*h)./(v+m2)).*(-h-k*v.*(v-a-1));
%normalizing vectors (to help plotting)
dv = dv./sqrt(dv.^2 + dh.^2);
dh = dh./sqrt(dv.^2 + dh.^2);
%generating the vector field
quiver(v, h, dv,dh,'AutoScaleFactor',0.5)
%changing axes limits, adding labels
axis([vmin vmax hmin hmax])
xlabel('$v$','Interpreter','latex')
ylabel('$h$','Interpreter','latex')

hold on
%plotting (v0,h0)=(0.5,0.2)
plot(0.5,0.2,'.black','MarkerSize', 20)
text(0.5,0.2,'(0.5,0.2)');
%solution to the system of diffeq's
f = @(t,x) [-k*x(1)*(x(1)-a)*(x(1)-1)-x(1)*x(2);(ep+((m1*x(2))/(m2+x(1)))*(-x(2)-8*x(1)*(x(1)-a-1)))];
[t,xa] = ode45(f,[0 500],[0.5 0.2]);
hold on
plot(xa(:,1),xa(:,2))
plot(xa(1,1),xa(1,2),'k*','LineWidth',2)
xlabel('v'), ylabel('h')

%figure3
close all;
clear all;
clc;
vmin = -0.05; vmax = 1.2;
hmin = -0.05; hmax = 2.75;
%step size for v and h
vstep = 0.1; hstep = 0.15;
%generating mesh for plotting
[v, h] = meshgrid(vmin:vstep:vmax, hmin:hstep:hmax);
%parameter values
a = 0.15; k = 8; ep = 0.002; m1 = 0.2; m2 = 0.3;
%system of equations dv/dt = f(v,h) and dh/dt = g(v,h)
dv = -k*v.*(v-a).*(v-1) - v.*h;
dh = (ep + m1*h./(v+m2)).*(-h-k*v.*(v-a-1));
%normalizing vectors (to help plotting)
dv = dv./sqrt(dv.^2 + dh.^2);
dh = dh./sqrt(dv.^2 + dh.^2);
%generating the vector field
quiver(v, h, dv,dh,'AutoScaleFactor',0.5)
%changing axes limits, add labels
axis([vmin vmax hmin hmax])
xlabel('$v$','Interpreter','latex')
ylabel('$h$','Interpreter','latex')
%plotting the solution to the system, with (v0,h0)=(0.5,0.2)
hold on
[v,h] = ode45(@system_ex,[-0.05,1.2],[0.5,0.2]);
plot(v(:,1),h(:,2))
%plotting (v0,h0)=(0.5,0.2)
plot(0.1,0.2,'.black','MarkerSize', 20)
text(0.1,0.2,'(0.5,0.2)');
%plotting the solution to system of diffeq's
function yprime = system_ex(t,y)
yprime = zeros(2,1);
yprime(1) = -8*y(1).*(y(1)-0.15).*(y(1)-1) - y(1).*y(2);
yprime(2) = (0.002 + 0.2*y(2)./(y(1)+0.3)).*(-y(2)-8*y(1).*(y(1)-0.15-1));
end

%figure4
k=8;
a=0.15;
mu1=0.2;
mu2=0.3;
sigma0=0.002;
tspan = [0 500];
f = @(t,x) [-k*x(1)*(x(1)-a)*(x(1)-1)-x(1)*x(2);
(sigma0+((mu1*x(2))/(mu2+x(1)))*(-x(2)-8*x(1)*(x(1)-a-1)))];
[t,xa] = ode45(f,tspan ,[0.25 0]);
hold on
plot(xa(:,1),xa(:,2))
v=-0.05:.01:1.2;
h=-0.05:.01:2.75;
ha = -k.*(v-a).*(v-1);
hb=((-sigma0 .*(v+mu2)) ./mu1);
hc= -k .*v .*(v-a-1);
x2=zeros(1,10);
y=linspace(-0.5,2.75,10);
axis = ([-0.05 1.2 -0.05 2.75]);
plot(v,ha, 'b');
hold on;
plot(v,hc,'r');
hold on;
plot(0,0,'k*')
plot(0,0,'ko')
hold on;
plot(v,hb,'r');
plot(x2,y,'b');
hold off;
legend('simulate solution','dv/dt','dh/dt','Equilibrium point');
xlabel('v')
ylabel('h')

%figure5
k=8;
a=0.15;
sigma0=0.002;
tspan = [0 500];
f = @(t,x) [-k*x(1)*(x(1)-a)*(x(1)-1)-x(1)*x(2);
(sigma0+((0.2*x(2))/(0.3+x(1)))*(-x(2)-8*x(1)*(x(1)-a-1)))];
[t,xa] = ode45(f,tspan ,[0.25 0]);
hold on
plot(t,xa(:,1))
plot(t,xa(:,2))
xlabel('t')
legend('v(t)','h(t)')

%figure6
function output = f(t)
T=100;
if ( mod(t,T) >= 10.0) && (mod(t,T) <= 13.0 )
output = 0.25;
else
output = 0;
end
end
%--------
function y = g(t,x)
mu1 = 0.2;
mu2=0.3;
sigma0=0.002;
v = x(1);
h = x(2);
y = [-8.*v.*(v-0.15).*(v-1)-v.*h+f(t);((sigma0+((mu1.*h)./(mu2+v))).*((-h)-8.*v.*(v-(0.15)-1)))];
end
%------
X0 = [0, 0];
tspan = 0:0.2:500;
[t,xa] = ode45(@g,tspan,X0);
plot(t,xa(:,1));
hold on
plot(t,xa(:,2));
xlabel('t')
legend('v(t)','h(t)')

%figure7
T=[100,90,80,70,60,50];
APD=[20.6,20.2,19.6,19,18,17];
scatter(T,APD);
xlabel('T')
ylabel('APD')
title('APD vs Stimulation Period T')

%figure8
T=[100,90,80,70,60,50];
h=[0.017442,0.020004,0.023303,0.027845,0.034229,0.043878];
scatter(T,h);
xlabel('T')
ylabel('steady-state h')
title('Steady-State h vs Stimulation Period T')

%figure9
APD=[20.6,20.2,19.6,19,18,17];
h=[0.017442,0.020004,0.023303,0.027845,0.034229,0.043878];
scatter(APD,h);
xlabel('APD')
ylabel('Steady-state h')
title('APD vs Steady-state h')