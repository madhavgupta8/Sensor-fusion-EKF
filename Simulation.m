clc; clear; 
close all;


% case 1: linear
n=8; %number of states
box_width = 2;
x0 =20; y0= -3; z0 = 1; T0=0; %initial position and orientation

a = -18.52; %deceleration m/s2

xdot0 = 22.5; %initial velocity %m/s
ydot0 = 0; %rate of change of height
zdot0 = 0;
Tdot0 = 0;

dt = 0.001;
t = 0:dt:1.5;
nt = length(t);

Xnonoise = [x0; xdot0; y0; ydot0; z0; zdot0; T0; Tdot0];

for k=1:(nt-1)
    xdot = Xnonoise(2,k);
    ydot = Xnonoise(4,k);
    zdot = Xnonoise(6,k);
    Tdot = Xnonoise(8,k)';
    Xnonoise(:, k+1)= Xnonoise(:,k) + dt*[xdot; a; ydot; 0; zdot; 0; Tdot; 0];
end

plot_birdseyeview(Xnonoise,[],[],'Truth: Birds Eye View', [1,3]);
PrepFigPresentation(gcf)
plot_birdseyeview ( Xnonoise ,[] ,[] , "True Trajectory XZ" ,[1,5])
PrepFigPresentation ( gcf )
camup ([0 1 0]) ;


% %case 2: rotating
% 
% n=8; %number of states
% box_width = 2;
% x0 = -8; y0= -3; z0 = 3.2; T0 = deg2rad(-10); %initial position and orientation
% 
% dt = 0.001;
% t = 0:dt:1;
% nt = length(t);
% 
% omega = -T0/t(end);
% omegadot = -omega/t(end);
% 
% a = -18.52; %deceleration m/s2
% vi = 22.5; %initial velocity m/s
% 
% xdot0 = vi*cos(T0); %initial velocity %m/s
% ydot0 = 0; %rate of change of height
% zdot0 = vi*sin(T0);
% Tdot0 = 0;
% 
% Xnonoise = [x0; xdot0; y0; ydot0; z0; zdot0; T0; Tdot0];
% 
% for k=1:(nt-1)
%     xdot = Xnonoise(2,k);
%     ydot = Xnonoise(4,k);
%     zdot = Xnonoise(6,k);
%     Tdot = Xnonoise(8,k);
%     T = Xnonoise(7,k);
%     Xnonoise(:, k+1)= Xnonoise(:,k) + dt*[xdot; a*cos(T); ydot; 0; zdot; a*sin(T); Tdot; omegadot];
% end
% 
% plot_birdseyeview(Xnonoise,[],[],'Truth: Birds Eye View', [1,3]);
% PrepFigPresentation(gcf)
% plot_birdseyeview(Xnonoise ,[] ,[] , "True Trajectory XZ",[1,5])
% PrepFigPresentation(gcf)
% camup ([0 1 0]) ;

function plot_birdseyeview (x1 ,x2 ,P2 , title_name , ii_plot );
% x1 is the true value or reference comparison
% x2 ,P2 is the estimator state and covariance
%
label_name ={ 'x','vx ','y','vy ','z','vz ','theta', 'thetadot'};
ii_x1 =[]; ii_x2 =[]; ii_P2 =[]; % for legend
figure ;
if ~ isempty (x1)
    plot (x1( ii_plot (1) ,:) ,x1( ii_plot (2) ,:) , 'Color' ,[0 0.5 0]) ; ii_x1 =1;
end
hold on;
if ~ isempty (x2),
    plot (x2( ii_plot (1) ,:) ,x2( ii_plot (2) ,:) ,'b-'); ii_x2 =2;
end
xlabel ( label_name ( ii_plot (1) )); ylabel ( label_name ( ii_plot (2) )); grid ;
hold off ;
legend_names ={ 'true trajectory','estimated trajectory','3\ sigma bound '};
legend ( legend_names { ii_x1 }, legend_names { ii_x2 }, legend_names { ii_P2 },'Location','South')
title (title_name);
PrepFigPresentation (gcf);
end

function PrepFigPresentation(fignum);
%
% prepares a figure for presentations
%
% Fontsize: 14
% Fontweight: bold
% LineWidth: 2
% 
figure(fignum);
fig_children=get(fignum,'children'); %find all sub-plots

for i=1:length(fig_children),
    
    set(fig_children(i),'FontSize',16);
    set(fig_children(i),'FontWeight','bold');
    
    fig_children_children=get(fig_children(i),'Children');
    set(fig_children_children,'LineWidth',2);
end
end