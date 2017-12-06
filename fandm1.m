%Calculation of forces and Moments

function out = fandm1(in)

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

%inputs
x = in(1);
x_dot = in(2);
y = in(3);
y_dot= in(4);
z = in(5);
z_dot= in(6);
roll = in(7);
roll_dot= in(8);
pitch = in(9);
pitch_dot= in(10);
yaw = in(11);
yaw_dot= in(12);
omega = in(13:16);


%Copied
V=sqrt(x_dot^2+y_dot^2);  % horizontal speed [m/s]
v1=sqrt(-0.5*V^2+sqrt((0.5*V^2)^2+(w/(2*rho*A))^2)); % Inflow velocity [m/s]
lambda=(v1+z_dot)/(OmegaH*R); % Inflow ratio [less]
mu=V/(OmegaH*R); % advance ratio [less]
muX=x_dot/(OmegaH*R); % advance ratio in x axis [less]
muY=y_dot/(OmegaH*R); % advance ratio in y axis [less]

%Coefficients
Ct=sigma_*a*(((1/6)+(mu^2)/4)*theta0-((1+mu^2)*thetatw/8)-lambda/4);
ChX=sigma_*a*((muX*Cd/(4*a))+(0.25*lambda*muX*(theta0-0.5*thetatw)));
ChY=sigma_*a*((muY*Cd/(4*a))+(0.25*lambda*muY*(theta0-0.5*thetatw)));
Cq=sigma_*a*((1/(8*a))*(1+mu^2)*Cd+lambda*((theta0/6)-(thetatw/8)-(lambda/4)));
CrX= -sigma_*a*(muX*(theta0/6-thetatw/8-lambda/8));
CrY= -sigma_*a*(muY*(theta0/6-thetatw/8-lambda/8));

T = zeros(4,1);
Hx = zeros(4,1);
Hy = zeros(4,1);
Q = zeros(4,1);
RRx = zeros(4,1);
RRy = zeros(4,1);

for(i=1:4)
T(i) = Ct*rho*A*((omega(i)*R)^2);           %Thrust Force
Hx(i) = ChX*rho*A*((omega(i)*R)^2);         %Hub force in x-direction
Hy(i) = ChY*rho*A*((omega(i)*R)^2);         %Hub force in y-direction 
Q(i) = Cq*rho*A*((omega(i)^2)*(R)^3);       %Drag  
RRx(i) = CrX*rho*A*((omega(i)^2)*(R)^3);    %Rolling moment x component
RRy(i) = CrY*rho*A*((omega(i)^2)*(R)^3);    %Rolling moment y component
end

%outputs
out(1:4)= T;
out(5:8)= Hx;
out(9:12)= Hy;
out(13:16)=Q;
out(17:20)= RRx;
out(21:24)= RRy;

