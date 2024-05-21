%% C FOR CAD %% MATLAB CODE for Coupler Curve Plotting  %%


%% Section 1 : Input Parameters
%--------------------------------------------------------------------------
%Link Lengths, Initial position angle and Coupler angle
%--------------------------------------------------------------------------
clc; clear; close all;
format SHORTG
a=20;                        %Crank 
b=100;                        %Coupler Link
c=50;                        %Output Link
d=100;                       %fixed Link
e=30;                        %CouplerExtension
alpha=90;                    %Coupler Extenstion Angle          
ar=alpha*pi/180;             %Degree to radian of Coupler Extenstion Angle
theta=10:1:270;              %Initial angle with limit of rotation 
T=theta*pi/60;              %Degree to radian of Initial angle
f=sqrt(a*a+d*d-2*a*d*cos(T)); 

%% Section 2 :Calculations 
%-------------------------- ------------------------------------------------
%Angles and Coupler Curve Coordinates 
%--------------------------------------------------------------------------
for i=1:length(theta) 
 %Calculation of Angles
G(i)=atan((a*sin(T(i)))/(d-a*cos(T(i))));
Bita(i)=acos((f(i).*f(1,i)+b*b-c*c)/(2*b*f(i)))-G(i);

%Coordinates 
X(i)=a*cos(T(i))+e*cos(ar+Bita(i));
Y(i)=a*sin(T(i))+e*sin(ar+Bita(i));


%% Section 3 :Plotting and Animation 
%-------------------------- ------------------------------------------------
%Plotting of Fourbar and Coupler Curve with trace 
hold on;
title('FourBar Coupler Curve Plotting')
xlim([-60,200])
ylim([-60,200])
plot(X,Y,'b','Linewidth',1)

Trace= viscircles([X(i) Y(i)],2,'Color','k');

Cx=[0 a*cos(T(i)) b*cos(Bita(i))+a*cos(T(i)) d a*cos(T(i))+e*cos(Bita(i)+ar)];
Cy=[0 a*sin(T(i)) b*sin(Bita(i))+a*sin(T(i)) 0 a*sin(T(i))+e*sin(Bita(i)+ar)];

crank=line([0 Cx(2)],[0 Cy(2)],'color','r','linewidth',5);%link 1
coupler=line([Cx(2) Cx(3)],[Cy(2) Cy(3)],'color','g','linewidth',5);%link 2
output=line([Cx(3) d],[Cy(3) 0],'color','y','linewidth',5);%link 3
fix=line([0 d],[0 0],'color','k','linewidth',5);%link4
eline=line([Cx(2) Cx(5)],[Cy(2) Cy(5)],'color','c','linewidth',3);
eoline=line([Cx(5) Cx(3)],[Cy(5) Cy(3)],'color','c','linewidth',3);%link e joint

%Animation of Plotting by Deleting the previous plotting 
pause(0.005)
delete(crank)
delete(coupler)
delete(output)
delete(eline)
delete(eoline)
delete(Trace)
end
%% Back to Initial Position 
crank=line([0 Cx(2)],[0 Cy(2)],'color','r','linewidth',5);%link 1
coupler=line([Cx(2) Cx(3)],[Cy(2) Cy(3)],'color','g','linewidth',5);%link 2
output=line([Cx(3) d],[Cy(3) 0],'color','y','linewidth',5);%link 3
fix=line([0 d],[0 0],'color','k','linewidth',5);%link4
eline=line([Cx(2) Cx(5)],[Cy(2) Cy(5)],'color','c','linewidth',3);
eoline=line([Cx(5) Cx(3)],[Cy(5) Cy(3)],'color','c','linewidth',3);%link e joint


% Velocity analysis
clc;
clear all;
close all;

format SHORTG
L2=20;                         
L3=50;                          
L4=100;                         
L1=100;                         
gama=180 ;                  
delta=(180-gama)/2;            
r=delta*pi/180;               
BE=2*L3*cos(r);                
theta=0:1:360;                
T=theta*pi/180;
 
k1=L1/L2;
k2=L1/L3;
k3=(L2^2 -L3^2 +L4^2 +L1^2)/(2*L2*L3);
k4=L1/L3;
k5=(L4^2 -L1^2 -L2^2 -L3^2)/(2*L2*L3);
 theta=0:1:360;              
T=theta*pi/180;
A= cos (T) -k1 - (k2*cos(T))+k3;
B=-2*sin(T);
C= k1 - ((k2+1)*cos(T)) + k3;
D=cos(T)-k1+k4*cos(T)+k5;
E=-2*sin(T);
F=k1+((k4-1)*cos(T))+k5;
x1=(-B-(sqrt((B.^2)-4*(A.*C)))).*(.5*(A.^-1));
x2=(-E-(sqrt((E.^2)-4*(D.*F)))).*(.5*(D.^-1));
omega2=(25*2*pi)/60;                           
 
for i=1:length(theta) ;   
theta2(i)=T(i);                  
    theta4(i)= 2*atan(x1(i));
     theta3(i)= 2*atan(x2(i));
     omega3(i)=((L2*omega2)/L3)*((sin(theta4(i)-theta2(i))/(sin(theta3(i)-theta4(i)))));
      omega4(i)=((L2*omega2)/L4)*((sin(theta2(i)-theta3(i))/(sin(theta4(i)-theta3(i)))));
       VA(i)=((((((L2*omega2)^2*((sin(theta2(i)))^2)+(j*cos(theta2(i)))^2)))));
      VPA(i)=(((BE*omega3(i))^2*((sin(theta3(i)+delta))^2)+(j*cos(theta3(i)+delta)))^2);
       VAx(i)=((L2*omega2)*(-sin(theta2(i))));
  VAy(i)=((((L2*omega2)*(cos(theta2(i))))));
  VPAx(i)=(((BE*omega3(i))*((-sin(theta3(i)+delta)))));
  VPAy(i)=((BE*omega3(i))*(cos(theta3(i)+delta)));
VPx(i)=VAx(i)+VPAx(i);
VPy(i)=VAy(i)+VPAy(i);
velocity(i)=sqrt((VPx(i)^2)+(VPy(i)^2));
 
hold on;
title('velocity VS. Theta 2 ')
figure(1);
plot(theta2,velocity,'color','g')
end

% Accelration analysis
clc; 
clear all; 
close all;
format SHORTG
L2=3;                         
L3=6.75;                          
L4=6.75;                         
L1=5.625;   
                      
BE=7.5;                       
delta=0;                      
r=delta*pi/180;               
theta=0:1:360;                 
T=theta*pi/180;               
 
k1=L1/L2;
k2=L1/L3;
k3=(L2^2 -L3^2 +L4^2 +L1^2)/(2*L2*L3);
k4=L1/L3;
k5=(L4^2 -L1^2 -L2^2 -L3^2)/(2*L2*L3);

A= cos (T) -k1 - (k2*cos(T))+k3;
B=-2*sin(T);
C= k1 - ((k2+1)*cos(T)) + k3;
D=cos(T)-k1+k4*cos(T)+k5;
E=-2*sin(T);
F=k1+((k4-1)*cos(T))+k5;
x1=(-B-(sqrt((B.^2)-4*(A.*C)))).*(.5*(A.^-1));
x2=(-E-(sqrt((E.^2)-4*(D.*F)))).*(.5*(D.^-1));

 omega2=(25*2*pi)/60;
 alpha2= 0;              
 
 
for i=1:length(theta) ;   
theta2(i)=T(i);                  
    theta4(i)= 2*atan(x1(i));
     theta3(i)= 2*atan(x2(i));
     omega3(i)=((L2*omega2)/L3)*((sin(theta4(i)-theta2(i))/(sin(theta3(i)-theta4(i)))));
      omega4(i)=((L2*omega2)/L4)*((sin(theta2(i)-theta3(i))/(sin(theta4(i)-theta3(i)))));
      VA(i)=((((((L2*omega2)^2*((-sin(theta2(i)))^2)+(j*cos(theta2(i)))^2)))));
      VPA(i)=(((BE*omega3(i))^2*((-sin(theta3(i)+delta))^2)+(j*cos(theta3(i)+delta)))^2);
 VAx(i)=((((((L2*omega2)*((-sin(theta2(i)))))))));
  VAy(i)=((((L2*omega2)*(cos(theta2(i))))));
  VPAx(i)=(((BE*omega3(i))*((-sin(theta3(i)+delta)))));
  VPAy(i)=((BE*omega3(i))*(cos(theta3(i)+delta)));
VPx(i)=VAx(i)+VPAx(i);
VPy(i)=VAy(i)+VPAy(i);



Aa=L4*sin(theta4(i));
Ba=L3*sin(theta3(i));
Ca=L2*alpha2*sin(theta2(i))+L2*(omega2^2)*cos(theta2(i))+L3*((omega3(i))^2)*cos(theta3(i))-L4*((omega4(i))^2)*cos(theta4(i));
Da=L4*cos(theta4(i));
Ea=L3*cos(theta3(i));
Fa=L2*alpha2*cos(theta2(i))-L2*((omega2)^2)*sin(theta2(i))-L3*((omega3(i))^2)*sin(theta3(i))+L4*((omega4(i))^2)*sin(theta4(i));
alpha3(i)=(Ca.*Da-Aa.*Fa)/(Aa.*Ea-Ba.*Da);
alpha4(i)=(Ca.*Ea-Ba.*Fa)/(Aa.*Ea-Ba.*Da);

AAx(i)=-L2*alpha2*sin(theta2(i))-L2*(omega2^2)*cos(theta2(i));
AAy(i)=L2*alpha2*cos(theta2(i))-L2*(omega2^2)*sin(theta2(i));
APAx(i)=(BE*alpha3(i)*-sin(theta3(i)+delta))-(BE*(omega3(i))^2*cos(theta3(i)+delta));
APAy(i)=(BE*alpha3(i)*cos(theta3(i)+delta))-(BE*(omega3(i))^2*sin(theta3(i)+delta));
APx(i)=APAx(i)+AAx(i);
APy(i)=APAy(i)+AAy(i);
Acceleration(i)=sqrt((APx(i)^2)+(APy(i)^2));

hold on;
title('Acceleration VS. Theta 2 ')  
figure(2);
   plot(theta2,Acceleration,'color','blue') ;
        
end


%% Force on each link analysis
% %% Mass of each link
% m=A*t*ro
%let width; thickness=0.5; Ro(steel)=2710*10^(-60);
width=2;thickness=0.5;
Ro=2710*10^(-6);
a=20;b=100;c=50;T4=0;delta=0
theta2=60;theta4=90;theta3=-42.31;
w2=1.047; w3=1.8250; w4=2.09439;
alpha2=0.5235;alpha3=0.90035;alpha4=1.7874;    
width=2;thickness=0.5;
Ro=2710*10^(-6);
m2=a*width*thickness*Ro;
m3=b*width*thickness*Ro;
m4=c*width*thickness*Ro;
%% Mass moment of inertia for each link
IG2=(1/12)*m2*(a^2+width^2);
IG3=(1/12)*m3*(4*b^2+width^2);
IG4=(1/12)*m4*(c^2+width^2);
%% Defining Parameters
R_CG2=a/2;
R_CG4=c/2;
R_CG3=b;
%% Force at point P
Fpx=20*cosd(45);
Fpy=20*sind(45);
T4=0;
%% X&Y Components of the position vectors
R_12x=R_CG2*cosd(theta2+180);
R_12y=R_CG2*sind(theta2+180);
R_32x=R_CG2*cosd(theta2);
R_32y=R_CG2*sind(theta2);
R_23x=R_CG3*cosd(theta3+180);
R_23y=R_CG3*sind(theta3+180);
R_43x=0;
R_43y=0;
R_34x=R_CG4*cosd(theta4);
R_34y=R_CG4*sind(theta4);
R_14x=R_CG4*cosd(theta4+180);
R_14y=R_CG4*sind(theta4+180);
R_px=b*cosd(theta3+delta);
R_py=b*sind(theta3+delta);
a_G2x=-R_CG2*alpha2*sind(theta2)-a*(w2^2)*cosd(theta2);
a_G2y=R_CG2*alpha2*cosd(theta2)-a*(w2^2)*sind(theta2);
a_G3x=-a*alpha2*sind(theta2)-a*(w2^2)*cosd(theta2)-alpha3*R_CG3*sind(theta3)-R_CG3*(w3^2)*cosd(theta3);
a_G3y=a*alpha2*cosd(theta2)-a*(w2^2)*sind(theta2)+alpha3*R_CG3*cosd(theta3)-R_CG3*(w3^2)*sind(theta3);
a_G4x=-R_CG4*alpha4*sind(theta4)-c*(w4^2)*cosd(theta4);
a_G4y=R_CG4*alpha4*cosd(theta4)-c*(w4^2)*sind(theta4);
% %% Solving 9x9 Matrix
AF=[1 0 1 0 0 0 0 0 0;
0 1 0 1 0 0 0 0 0;
-R_12y R_12x -R_32y R_32x, 0 0 0 0 1;
0 0 -1 0 1 0 0 0 0;
0 0 0 -1 0 1 0 0 0;
0 0 R_23y -R_23x -R_43y ,R_43x 0 0 0;
0 0 0 0 -1 0 1 0 0;
0 0 0 0 0 -1 0 1 0;
0 0 0 0 R_34y -R_43x -R_14y R_14x 0];
clc
Z=[m2*a_G2x; m2*a_G2y; IG2*alpha2; m3*a_G3x-Fpx; m3*a_G3y-Fpy;IG3*alpha3-R_px*Fpy+R_py*Fpx; m4*a_G4x; m4*a_G4y; IG4*alpha4-T4];
Fo=AF\Z;
Fo=10^(-2)*Fo; %to make it in N
F12x = Fo(1);
F12y = Fo(2);
F32x = Fo(3);
F32y = Fo(4);
F43x = Fo(5);
F43y = Fo(6);
F14x = Fo(7);
F14y = Fo(8);
T12 = Fo(9)








