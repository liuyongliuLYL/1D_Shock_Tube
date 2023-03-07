clc;
%clear;
u=load("velocity.dat");
rou=load("desnity.dat");
p=load("pressure.dat");
% c=load("c.dat");
% ave_rou=load("ave_rou.dat");
% ave_p=load("ave_p.dat");
% ave_u=load("ave_u.dat");
% ave_c=load("ave_c.dat");
%u_exact=load("../Sod_exact/velocity.dat");
x=linspace(-0.5,0.5,101);

figure(1);
ceng=130;
%plot(xSpn, pSpn, xSpn, uSpn, xSpn, rSpn);
hold on
plot(x,u(:,ceng),x,rou(:,ceng),x,p(:,ceng));

% xlim([0 1]);
% ylim([-0.2 1.5]);
%title('t=0.30时刻');


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
% loops = 401;
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