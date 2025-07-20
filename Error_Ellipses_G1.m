clc
clear all
close all
format long g

disp('Error Ellipses Project')

syms x2 y2 x3 y3 x4 y4 x5 y5 x6 y6
%% Initial amounts...
M = load('Observation.txt');
x1=1000;
y1=1000;
L = M(:,11); % Average Observations
L(13,1) = deg2rad(96.3175);
%  Initial coordinates amounts
X0=[1361.253,960.006,1550.110,1099.304,1297.240,1484.943,1093.298,1416.196,1191.891,1159.498];
% Dis equations
X=[x2, y2, x3, y3, x4 ,y4, x5, y5, x6, y6];
  F(1)  = sqrt((x2-x1)^2+(y2-y1)^2);
  F(2)  = sqrt((x3-x1)^2+(y3-y1)^2);
  F(3)  = sqrt((x5-x1)^2+(y5-y1)^2);
  F(4)  = sqrt((x3-x2)^2+(y3-y2)^2);
  F(5)  = sqrt((x4-x2)^2+(y4-y2)^2); 
  F(6)  = sqrt((x6-x2)^2+(y6-y2)^2);
  F(7)  = sqrt((x4-x3)^2+(y4-y3)^2);
  F(8)  = sqrt((x6-x3)^2+(y6-y3)^2);
  F(9)  = sqrt((x4-x6)^2+(y4-y6)^2);
  F(10) = sqrt((x5-x6)^2+(y5-y6)^2);
  F(11) = sqrt((x5-x4)^2+(y5-y4)^2);
  F(12) = sqrt((x3-x5)^2+(y3-y5)^2);
  F(13) = pi-atan(abs((x2-x1)/(y2-y1)));
  F=F';
  
%% Observations Var-Cov matrix 

% by camera constant error
Qy1 = zeros(13); % faster...
for i=1:12
    Qy1(i,i) = ((1.5)^2 + (2*L(i)/1000)^2)*10^-7;
end
Qy1(13,13) = (deg2rad(1/3600))^2;
Qy1 = inv(Qy1);

% by variance calculation
Qy2 = zeros(13); % faster...
for i=1:12
    Qy2(i,i) = var(M(i,1:10))/10;
end
Qy2(13,13) = (deg2rad(1/3600))^2;
Qy2 =inv(Qy2); 

%% Proration

A = jacobian(F,X);     
A0 = eval(subs(A,X,X0));   
L0 = eval(subs(F,X,X0));   
dL = L - L0;
dx = inv((A0)'*Qy2*(A0))*(A0)'*Qy2*dL;
dX00 = (dx)';
counter=1;
while abs(dx)>10^-6
    X0     = X0+dX00;
    A0       = eval(subs(A,X,X0));
    L0     = eval(subs(F,X,X0));
    dL     = L - L0;
    dx     = inv((A0)'*Qy2*(A0))*(A0)'*Qy2*dL;
    dX00 = (dx)';
    counter=counter+1;
end
    X0     = X0+dX00;
    
%% Remains

V_prim=A0*dx-dL; 

%% Degrees of freedom

dr=length(M)-length(X);
answer1=[' Degrees of freedom = ', num2str(dr)];
disp(answer1)

%% 2end Variance Factor 

S=(V_prim'*Qy2*V_prim)/dr; 
answer2=['   2end Variance Factor = ', num2str(S)];
disp(answer2)

Limitation1=chi2inv(0.025,dr);
Limitation2=chi2inv(0.975,dr);

%% Precision of parameters

Qx = inv(A0'*Qy2*A0); 

%% test chi-square

if S*dr>Limitation1 && S*dr<Limitation2
    disp('    Test Result :  Accepted ');
else
    disp('    Test Result :  Rejected! ')
    %% Estimation of Precision of parameters
    
    Qx=Qx*S;
    
end

%% Stations & Vectors plot

% Stations
Vx(1,1)=1000; %first_x
Vy(1,1)=1000; %first_y
for i = 1:5
    Vx(i+1,1) = X0(1,2*i-1);
    Vy(i+1,1) = X0(1,2*i);
end

figure(1)
plot(Vx(1,1),Vy(1,1),'^','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b')
hold on
plot(Vx(2:6,1),Vy(2:6,1),'s','LineWidth',3,'MarkerSize',8,'MarkerEdgeColor','r');
title('Error Ellipses','fontsize',20)
hold on 

% Vectors
for i=1:12
    s1=M(i,12);
    s2=M(i,13);
    Vx_1=[Vx(s1,1),Vx(s2,1)];
    Vy_1=[Vy(s1,1),Vy(s2,1)];
    line(Vx_1,Vy_1);
end

axis equal
 hold on
 
 %% Absolute Precision of points
 
 semi_majorAxis=zeros(1,6);
 semi_minorAxis=zeros(1,6);
 t_exact=zeros(1,6);
 Alpha_exact=zeros(1,6);

for i=1:5
    QxA=Qx(2*i-1:2*i,2*i-1:2*i);
    semi_majorAxis(i+1)= sqrt(max(eig(QxA)));
    semi_minorAxis(i+1)= sqrt(min(eig(QxA)));
    t_exact(i+1)=(atan(2*QxA(1,2)/(QxA(1,1)-QxA(2,2))));
    if t_exact(i+1)<0
        if QxA(1,2)>0
            t_exact(i+1)= t_exact(i+1)+180;%2
        elseif QxA(1,2)<0
            t_exact(i+1)= t_exact(i+1)+360;%4
        end
    end
    if t_exact(i+1)>0
        if QxA(1,2)<0
            t_exact(i+1)= t_exact(i+1)+180;%3
        elseif QxA(1,2)>0
            t_exact(i+1)= t_exact(i+1);%1
        end
    end
    Alpha_exact(i+1)=t_exact(i+1)/2;
    ellipse(semi_majorAxis(i+1)*20000,semi_minorAxis(i+1)*20000,Alpha_exact(i+1),Vx(i+1,1),Vy(i+1,1),'g');
    axis equal
    hold on
end


%% Comparitive Precision of Points
alpha_com=zeros(1,12);
for i=1:12
    s1=M(i,12);
    s2=M(i,13);
    
    if s1==1    % Exception!...
        continue
    end
    
    % 
    QxC(1:2,1:2)=Qx(2*s1-3:2*s1-2,2*s1-3:2*s1-2);
    QxC(3:4,3:4)=Qx(2*s2-3:2*s2-2,2*s2-3:2*s2-2);
    QxC(1:2,3:4)=Qx(2*s1-3:2*s1-2,2*s2-3:2*s2-2);
    QxC(3:4,1:2)=QxC(1:2,3:4)'; % error c_v
    A=[-1 0 1 0;0 -1 0 1];
    Qxn=A*QxC*A'; % error affect between two points
    
    a_com(i)= sqrt(max(eig(Qxn)));
    b_com(i)= sqrt(min(eig(Qxn)));
    t_com(i)=(atan(2*Qxn(1,2)/(Qxn(1,1)-Qxn(2,2))));
    if  t_com(i)<0
        if Qxn(1,2)>0
            t_com(i)= t_com(i)+180;
        elseif Qxn(1,2)<0
            t_com(i)= t_com(i)+360;
        end
    end
    if  t_com(i)>0
        if Qxn(1,2)<0
            t_com(i)= t_com(i)+180;
        end
    end
    
    alpha_com(i)=t_com(i)/2;
    x_ave=(Vx(s1,1)+Vx(s2,1))/2;
    y_ave=(Vy(s1,1)+Vy(s2,1))/2;
    plot(x_ave,y_ave,'d');
    ellipse(20000*a_com(i),20000*b_com(i),alpha_com(i),x_ave,y_ave,'m')
    axis equal
    hold on
    
end

Absolute_surface=pi*semi_majorAxis.*semi_minorAxis;
compaitive_surface=pi*a_com.*b_com;

%% Surface values comprarison

figure(2)
plot(compaitive_surface(1,3:12),'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b')
hold on
plot(Absolute_surface,'r')
title('Absolute & Comparitive surface change','fontsize',15)
legend('compaitive Surface','Absolute Surface')

