clear;
%clc;
u=load("velocity.dat");
rou=load("desnity.dat");
p=load("pressure.dat");
x=linspace(-2,2,401);

figure(1);
plot(x,u(:,50),x,rou(:,50),x,p(:,50));
legend('速度','密度','压力');
xlim([-2 2]);
ylim([-0.2 1.2]);
title('t=0.5时刻');



u1=u(:,1);
rou1=rou(:,1);
p1=p(:,1);

h = figure;
Z = plot(x,u1,x,rou1,x,p1); legend('速度','密度','压力');
% axis tight manual%设置坐标轴范围和纵横比
% ax = gca;
% ax.NextPlot = 'replaceChildren';

loops = 101;
M(loops) = struct('cdata',[],'colormap',[]);

h.Visible = 'off';
for j = 1:loops
    u1=u(:,j);
    rou1=rou(:,j);
    p1=p(:,j);
    plot(x,u1,x,rou1,x,p1); legend('速度','密度','压力');
    drawnow
    M(j) = getframe;
end

h.Visible = 'on';

movie(M);