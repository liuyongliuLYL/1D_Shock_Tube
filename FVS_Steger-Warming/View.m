clear;
clc;
u=load("velocity.dat");
rou=load("desnity.dat");
p=load("pressure.dat");
x=linspace(-0.5,0.5,501);

figure(1);
ceng=2;
% scatter(x,u(:,ceng));
% figure(2);
% scatter(x,rou(:,ceng));
% figure(3);
% scatter(x,p(:,ceng));
%plot(xSpn, pSpn, xSpn, uSpn, xSpn, rSpn);hold on
plot(x,u(:,ceng),x,rou(:,ceng),x,p(:,ceng)); 
legend("速度","密度","压力");
% xlim([0 1]);
% ylim([-0.2 1.2]);
%title('t=0.30时刻');
% 
% 
% u1=u(:,1);
% rou1=rou(:,1);
% p1=p(:,1);
% 
% h = figure(2);
% Z = plot(x,u1,x,rou1,x,p1); legend('速度','密度','压力');
% % axis tight manual%设置坐标轴范围和纵横比
% % ax = gca;
% % ax.NextPlot = 'replaceChildren';
% 
% loops = 2001;
% M(loops) = struct('cdata',[],'colormap',[]);
% 
% h.Visible = 'off';
% for j = 1:loops
%     u1=u(:,j);
%     rou1=rou(:,j);
%     p1=p(:,j);
%     plot(x,u1,x,rou1,x,p1); legend('速度','密度','压力');
%     drawnow
%     M(j) = getframe;
% end
% 
% h.Visible = 'on';
% 
% movie(M);