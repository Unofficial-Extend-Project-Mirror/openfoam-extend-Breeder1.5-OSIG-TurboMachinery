%
% Radial profile at the inlet (-500mm) based on the Gyllenram and Nilsson
% velocity profile (Yokohama 2006) and on the k and epsilon Case1 default
% values
%
% Radial velocity = 0 [m/s]
% Turbulence Kinetic Energy = 2.0184 [m^2/s^2]
% Turbulence Eddy Dissipation = k^1.5 / l_t = 896.10908 [m^2/s^3]
%
% This case is part of the OpenFOAM TurboMachinery Special Interest Group
%
% For more information about the TurboMachinery SIG, please visit
% http://openfoamwiki.net/index.php/Sig_Turbomachinery
%
% 2008-04-10 Omar Bounous, Chalmers University of Technology
%

close all;
clear all;

N=30;
r=linspace (0,0.13,N);      % r-coordinate from the centerline
R=0.13;                     % Honeycomb radius
h=R/(N-1);

U_b=11.6;                     % Bulk velocity
V_z0r0_on_U_b=0.98;
n=4;

C=((R^2)*(1-V_z0r0_on_U_b)*(n+2))/(2*R^(n+2));

for i=1:N
V_ax(i)=U_b*(V_z0r0_on_U_b+C*(r(i))^n);
V_rad(i)=0;
V_circ(i)=0.59*U_b*((r(i)/0.13));
p(i)=0;
k(i)=2.0184;
epsilon(i)=896.10908;
end

Inlet_BC(N,7)=zeros;
Inlet_BC= [r' V_ax' V_rad' V_circ' p' k' epsilon'];

fid = fopen('Inlet_BC.csv','wt');
fprintf(fid,'#\n');
fprintf(fid,'# Radial profile at the inlet (-500mm) based on the Gyllenram and Nilsson velocity profile (Yokohama 2006) and on the k and epsilon Case1 default values\n');
fprintf(fid,'#\n');
fprintf(fid,'# Radial velocity = 0 [m/s]\n');
fprintf(fid,'# Turbulence Kinetic Energy = 2.0184 [m^2/s^2]\n');
fprintf(fid,'# Turbulence Eddy Dissipation = k^1.5 / l_t = 896.10908 [m^2/s^3]\n');
fprintf(fid,'#\n');
fprintf(fid,'# This case is part of the OpenFOAM TurboMachinery Special Interest Group\n');
fprintf(fid,'#\n');
fprintf(fid,'# For more information about the TurboMachinery SIG, please visit\n');
fprintf(fid,'# http://openfoamwiki.net/index.php/Sig_Turbomachinery\n');
fprintf(fid,'#\n');
fprintf(fid,'# 2008-04-10 Omar Bounous, Chalmers University of Technology\n');
fprintf(fid,'#\n');
fprintf(fid,'[Data]\n');
fprintf(fid,' R [ m ], Velocity Axial [ m s^-1 ], Velocity Radial [ m s^-1 ], Velocity Circumferential [ m s^-1 ], Pressure [ Pa ], Turbulence Kinetic Energy [ m^2 s^-2 ], Turbulence Eddy Dissipation [ m^2 s^-3 ]\n');
fprintf(fid,'%12.8f, %12.8f, %12.0f, %12.4f, %12.0f, %12.5f, %12.6f\n',Inlet_BC');
fclose(fid);

figure(1)
plot(r,V_ax)

figure(2)
plot(r,V_circ)
