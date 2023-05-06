%% MAE 6760 Final Project

clc; clear;
close all;

%% Simulation
n=8; %number of states
box_width = 2;
x0 = -20; y0= -3; z0 = 1; T0=0; %initial position and orientation

a = -12.5; %deceleration m/s2

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
    Xnonoise(:, k+1)= Xnonoise(:,k) + dt*...
        [xdot; 
        a; 
        ydot; 
        0; 
        zdot; 
        0; 
        Tdot; 
        0];
end

plot_birdseyeview(Xnonoise,[],[],'Truth: Birds Eye View', [1,3]);
PrepFigPresentation(gcf)
plot_birdseyeview ( Xnonoise ,[] ,[] , "True Trajectory XZ" ,[1,5])
PrepFigPresentation(gcf)
camup([0 1 0]) ;

%% cameras setup
psi1 = atan(Xnonoise(5 ,:)./Xnonoise(1 ,:));
phi1 = atan(Xnonoise(3 ,:)./sqrt(Xnonoise(1 ,:).^2 + Xnonoise(5 ,:).^2) );
psi2 = atan((box_width - Xnonoise (5 ,:))./ Xnonoise (1 ,:));
phi2 = atan(Xnonoise (3 ,:)./sqrt(Xnonoise(1,:).^2 + (box_width-Xnonoise(5 ,:)).^2));
% psi1 = atan2(Xnonoise(5 ,:), Xnonoise(1 ,:));
% phi1 = atan2(Xnonoise(3 ,:), sqrt(Xnonoise(1 ,:).^2 + Xnonoise(5 ,:).^2) );
% psi2 = atan2((box_width - Xnonoise (5 ,:)), Xnonoise (1 ,:));
% phi2 = atan2(Xnonoise (3 ,:), sqrt(Xnonoise(1,:).^2 + (box_width-Xnonoise(5 ,:)).^2));
nz = 4;
nw = 4;
R= eye(nz)*0.005^2; v= sqrtm(R)*randn(nz,nt);
Q=eye(nw)*0.05^2; w= sqrtm(Q)*randn (nw ,nt);
Z = (180/pi)*[psi1; phi1; psi2; phi2] +v;
Zacc = w;
figure(7)
plot(t, psi1)
hold on
plot(t, phi1)
plot(t, psi2, "g--")
plot(t, phi2, "--")
legend(["Psi1", "Phi1", "Psi2", "Phi2"])
ylabel("angle (radians)")
xlabel("time")
title("Measurements verses time")



%% EKF - linear case
x0 =[-20; 0; -3; 0; 1; 0; 0; 0];
P0=diag([1 1 0.001 0.001 1 1 0.001 0.001])*1^2;
%P0 = eye(n)*1^2;
n= length(x0);
xhatp =x0;Pp (1:n ,1:n ,1)=P0;
xhatu =x0;Pu (1:n ,1:n ,1)=P0;
h = 3; %height of camera

inn_T = zeros(nz,1);

for k=1:(nt-1)
    %predict state
    xhatp(:,k+1) = predict_state(xhatu(:,k), Zacc(:,k), dt, a);
    %predict covariance
    [F,G] = getFG(xhatu(:,k), dt);
    Pp(1:n, 1:n, k+1) = F*Pu(1:n, 1:n, k)*F' + G*Q*G';
    %Kalman Gain
    H = getH(xhatp(:, k+1));
    Zp = (180/pi)*...
        [atan(xhatp(5,k+1)/xhatp(1,k+1));
        atan(xhatp(3,k+1)/sqrt(xhatp(1,k+1)^2 + xhatp(5,k+1)^2));
        atan((box_width - xhatp(5,k+1))/xhatp(1,k+1));
        atan(xhatp(3,k+1)/sqrt(xhatp(1,k+1)^2 + (box_width - xhatp(5,k+1))^2))]; 
    inn = Z(:, k+1) - Zp;
    inn_T(:, k+1) = inn;
    S = H*Pp(1:n, 1:n, k+1)*H' + R;
    K = Pp(1:n, 1:n, k+1)*H'*inv(S);

    %Update
    xhatu(:, k+1) = xhatp(:, k+1) + K*(Z(:,k+1)-Zp);
    Pu(1:n, 1:n, k+1) = (eye(n) - K*H)*Pp(1:n, 1:n, k+1)*(eye(n) - K*H)' + K*R*K';

    %test hypothesis
    goal = check(xhatu(:, k+1), Pu(:, :, k+1));
end


plot_birdseyeview(Xnonoise, xhatu, Pu,'Estimated Trajectory XY', [1 ,3])
PrepFigPresentation(gcf)
plot_birdseyeview(Xnonoise, xhatu, Pu,'Estimated Trajectory XZ', [1 ,5])
PrepFigPresentation(gcf)
ii_plot = [1 5];
plot_estimator_error(t, Xnonoise, xhatu, Pu, ii_plot, "Estimator error : X and Y axis");

figure(6)
plot(t, inn_T(1,:))
title('Innovation(Z - Z_p) vs time: for measurement 1 (pan and rot)')

% %% Plotting Trajectory
% figure(6)
% for k =1:(nt)
%     goal = check(xhatu(:,k),Pu(:,:,k));
%     scatter3(Xnonoise(1,k),Xnonoise(3,k), Xnonoise(5,k),'green')
%     if goal == 1
%         scatter3(xhatu(1,k),xhatu(3,k),xhatu(5,k),'magenta')
%         hold on;
%     else
%         scatter3(xhatu(1,k),xhatu(3,k),xhatu(5,k),'blue')
%         hold on;
%     end
% end
% title('Trajectory of the car')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% camup([0 1 0])

%% EKF function calls

function Xkp1 = predict_state(Xk, U, dt, a)
wx = U(1); wy = U(2); wz = U(3); wT = U(4);
Xkp1 = Xk + dt*...
    [Xk(2);
    wx + a;
    Xk(4);
    wy;
    Xk(6);
    wz;
    Xk(8);
    wT];
end

function [F,G] = getFG(X,dt)
F = [1 dt 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 1 dt 0 0;
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 1 dt;
    0 0 0 0 0 0 0 1];
G = [0 0 0 0;
    dt 0 0 0;
    0 0 0 0;
    0 dt 0 0;
    0 0 0 0;
    0 0 dt 0;
    0 0 0 0;
    0 0 0 dt];
end

function H = getH(X)
L = 2; %2 m wide pit box
H = 5; %camera is placed 5 m above
x = X(1); y = X(3); z=X(5);T = X(7);

H = (180/pi)*...
    [-z/(x^2+z^2), 0, 0, 0, x/(x^2+z^2), 0, 0, 0;
    -x*y/((x^2+ y^2+ z^2)*sqrt(x^2+ z^2)), 0, sqrt(x^2+z^2)/(x^2+ y^2+ z^2), 0, -z*y/((x^2+y^2+z^2)*sqrt(x^2+z^2)), 0, 0, 0;
    -(L-z)/(x^2+(L-z)^2), 0, 0, 0, -x/(x^2+(L-z)^2), 0, 0, 0;
    -x*y/((x^2+ y^2+(L-z)^2)*sqrt(x^2+(L-z)^2)), 0, sqrt(x^2+(L-z)^2)/(x^2+ y^2+(L-z)^2), 0, 0, (L-z)*y/((x^2+ y^2+(L-z)^2)*sqrt(x^2+(L-z)^2)), 0, 0];
end



%% Hypothesis Testing Function

function goal = check(Xk, Pk)
prob = 1 - cdf('Normal', 1, Xk(5), sqrt(Pk(3,3)));
if prob>0.7
    goal = 1;
else
    goal = 0;
end
end

%% Plotting Functions

function plot_birdseyeview (x1 ,x2 ,P2 , title_name , ii_plot);
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

function [Xe,Ye] = calculateEllipseCov(X, P, nsig, steps) 
    %# This functions returns points to draw an ellipse 
    %# 
    %#  @param X     x,y coordinates 
    %#  @param P     covariance matrix 
    %# 
 
    error(nargchk(2, 4, nargin)); 
    if nargin<3, nsig = 1; end 
    if nargin<4, steps = 36; end 
    
    [U,S,V]=svd(P);
    s1=sqrt(S(1,1));s2=sqrt(S(2,2));angle=acos(U(1,1))*180/pi;
    x=X(1);
    y=X(2);

    %scale by nsig
    s1=nsig*s1;
    s2=nsig*s2;

    beta = angle * (pi / 180); 
    sinbeta = sin(beta); 
    cosbeta = cos(beta); 
 
    alpha = linspace(0, 360, steps)' .* (pi / 180); 
    sinalpha = sin(alpha); 
    cosalpha = cos(alpha); 
 
    Xe = x + (s1 * cosalpha * cosbeta - s2 * sinalpha * sinbeta); 
    Ye = y + (s1 * cosalpha * sinbeta + s2 * sinalpha * cosbeta); 
 
end

function plot_estimator(k,x1,x2,z,P2,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% z is the measurement
%
% DOES NOT CHECK FOR SIZES!!
%
figure;
if ~isempty(z),
    pp=plot(k,x1,'r-',k,x2,'b-.',k,z,'g:');set(pp(3),'color',[0 0.5 0]);
else,
    pp=plot(k,x1,'r-',k,x2,'b-.');
end
hold on;
if ~isempty(P2)
    plot(k,x2-2*sqrt(P2),'b:',k,x2+2*sqrt(P2),'b:');
    axis_name='state estimate';
    legend('true state','estimate','measurement','2\sigma bound','Location','Northeast');
else,
    axis_name='state';    
    legend('true state','no noise','measurement','Location','Northeast');    
end
hold off
xlabel('time (sec)');ylabel(axis_name);grid;
title(title_name);
PrepfigPresentation(gcf);
end

function plot_estimator_error(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'x','xdot','y','ydot','z','zdot', 'T', 'Tdot'};
figure;subplot(122);
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    err=x2(ii,:)-x1(ii,:);
    plot(t,err,'b-');
    hold on;
    if ~isempty(P2)
        tbound=[t fliplr(t)];
        xbound=[[err+2*sqrt(squeeze(P2(ii,ii,:)))'] fliplr([err-2*sqrt(squeeze(P2(ii,ii,:)))'])];
        patch(tbound,xbound,'b','EdgeColor','b','FaceAlpha',0.2,'EdgeAlpha',0.2);
        plot(t,zeros(length(t),1),'r--');
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 1.5]);set(gca,'xtick',[0:0.1:1.5]);
end
sgtitle(title_name);
%
PrepFigPresentation(gcf);
legend('estimator error','2\sigma bound','zero error','Location','South','fontsize',12);
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
