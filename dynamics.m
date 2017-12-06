% dynamics

function out = dynamics(in)

%///////////////////

% constants file
% INPUT: NaF
% OUTPUT: all constants

sp= 0.01; % [s] sampling period only used to compute YiA term in dinamica.m
g= 9.806;    % [m/s^2]
b= 3.13E-5;       % [N.s2] thrust factor in hover
d= 7.5E-7;         % [Nm.s2] drag factor in hover

%***** DEG2RAD ***** 
deg2rad=1.745329E-2; % [rad/deg]
rad2deg=57.295779;  % [deg/rad]

%***** Number of propellers  *****/

P=4;
%***** Arm length *****/
L= 0.232; 			%// [m]

%***** Total mass *****/
m = 0.53;	%// [kg]

%***** Inertia components *****/
Ixx = 6.228E-3;	%// [kg.m2]    
Iyy = 6.228E-3;	%// [kg.m2]
Izz = 1.121E-2;	%// [kg.m2]

%***** Rotor inertia *****/
r=4;    % reduction ratio
jm= 4e-7;	%// [kg.m2];  
jp= 6e-5;	%// [kg.m2];  
jr = jp+jm/r;	%// [kg.m2];

%***** CoG position *****/
% h= 0.058; 			%// [m] VERTICAL DISTANCE between CoG and propellers plan
h= 0; 			%// [m] more stable with h=0 !!!

%***** constants of the linear curve Omega=f(bin) *****
slo=2.7542; % slope (of the linear curve om=f(bin))
shi=3.627;  % shift

% *************************** AERODYNAMICS ***************************

% ********* Propeller ********* 
N=2; % number of blades
R=0.15; % [m] propeller radius
A=pi*(R^2); % [m^2] prop disk area
c=0.0394; % [m]   chord
theta0=0.2618;   % [rad] pitch of incidence
thetatw=0.045;   % [rad] twist pitch
sigma_=(N*c)/(pi*R); % solidity ratio (rotor fill ratio) [rad^-1]
a=5.7;  % Lift slope (given in literature)
Cd=0.052;   % Airfoil drag coefficient -- found by tests
Ac=0.005;  % [m^2] helicopter center hub area 

% ********* Gaz ********* 
rho= 1.293;     % [kg/m^3]   air density
nu=1.8e-5;      % [Pa.s] air viscosity @ 20deg

% ********* misc ********* 
w=(m/P)*g;   % [N] weigth of the helicopter/number of propellers
OmegaH=sqrt(w/b); % [rad/s] prop spd at hover
OmegaMax = 600; %[rad/s]
% ********* Longitudinal Drag Coefficients *********
Cx=1.32;    % estimation taken from literature [less]
Cy=1.32;    % estimation taken from literature [less]
Cz=1.32;    % estimation taken from literature [less]

% ***** OS4 volume *****
Vol=3.04E-4; 	% [m3]
PArchim = rho*g*Vol;  % [N]


% ***** ANNEXE ! ******
%Ftb=0.5*Cz*A*rho*v^2   % turbulent
%Fl=16*0.3*nu*v   % laminary 





%////////////


%input states and omegas
x = in(1);
x_dot = in(2);
y = in(3);
y_dot = in(4);
z = in(5);
z_dot= in(6);
roll = in(7);
roll_dot= in(8);
pitch = in(9);
pitch_dot= in(10);
yaw = in(11);
yaw_dot= in(12);
omega = in(13:16);

%input forces and moments
T=in(17:20);             %Thrusts [N]
Hx=in(21:24);            %Hub forces x component
Hy=in(25:28);            %HUb forces y component 
Q=in(29:32);             %Drag
RRx=in(33:36);           %Rolling moment x component 
RRy=in(37:40);           %Rolling moment x component

Om = omega(1) - omega(2) + omega(3) - omega(4);
global Om_old;

% Rolling moments
Rbg = -pitch_dot*yaw_dot*(Iyy-Izz);                
Rgp = jr*pitch_dot*Om;                            
Raa = L*(-T(2)+T(4));                           
Rhm = (Hy(1)+Hy(2)+Hy(3)+Hy(4))*h;              
Rrm = +RRx(1)-RRx(2)+RRx(3)-RRx(4);             
Rfm = 0.5*Cz*A*rho*roll_dot*abs(roll_dot)*L*(P/2)*L;   

% Pitch moments
Pgb = roll_dot*yaw_dot*(Izz-Ixx); 
Pgp = jr*roll_dot*Om; 
Paa = L*(-T(1)+T(3)); 
Phf = (Hx(1)+Hx(2)+Hx(3)+Hx(4))*h; 
Prm = +RRy(1)-RRy(2)+RRy(3)-RRy(4);            
Pfm = 0.5*Cz*A*rho*pitch_dot*abs(pitch_dot)*L*(P/2)*L; 

% Yaw moments
Ygb = pitch_dot*roll_dot*(Ixx-Iyy); 
Yict = jr*(Om-Om_old)/sp;           
Yct = +Q(1)-Q(2)+Q(3)-Q(4);              
Yhfx = (-Hx(2)+Hx(4))*L;                    
Yhfy = (-Hy(1)+Hy(3))*L; 

Om_old=Om;


% X forces
Xaa = (sin(yaw)*sin(roll)+cos(yaw)*sin(pitch)*cos(roll))*(T(1)+T(2)+T(3)+T(4)); 
Xdf = 0.5*Cx*Ac*rho*x_dot*abs(x_dot); 
Xhf = Hx(1)+Hx(2)+Hx(3)+Hx(4);

 % Y forces
Yaa = (-cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll))*(T(1)+T(2)+T(3)+T(4)); 
Ydf = 0.5*Cy*Ac*rho*y_dot*abs(y_dot); 
Yhf = Hy(1)+Hy(2)+Hy(3)+Hy(4); 

% Z forces
Zaa = (cos(pitch)*cos(roll))*(T(1)+T(2)+T(3)+T(4));              
Zaf = 0.5*Cz*A*rho*z_dot*abs(z_dot)*P + 0.5*Cz*Ac*rho*z_dot*abs(z_dot); 

%outputs
out(1) = roll_dot;
out(2) = (Rbg + Rgp + Raa + Rhm + Rrm - Rfm) /Ixx;
out(3) = pitch_dot;
out(4) = (Pgb - Pgp + Paa - Phf + Prm - Pfm) /Iyy;
out(5) = yaw_dot;
out(6) = (Ygb + Yict +Yct + Yhfx + Yhfy) /Izz;
out(7) = x_dot;
out(8) = (Xaa - Xdf - Xhf)/m;
out(9) = y_dot;
out(10) = (Yaa - Ydf - Yhf) / m;
out(11) = z_dot;
out(12) = -g + (Zaa - Zaf)/m;

