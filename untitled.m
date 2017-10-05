clc
clear all
close all



% Input Data %
fin_vel=csvread('FINAL_VEL'); %% read input file
in_vel=csvread('IN_VEL'); %% read input file
En=dlmread('ENERGY'); %% read input file

% Storing particles' energies %
for i=1:1:(length(En))   %% Loop to store Potential, Kinetic and Total energies.
    E_Pot(i,:)=En(i,3);
    E_Kin(i,:)=En(i,2);             %% also to store time stamps
    E_Tot(i,:)=En(i,4);
    time_E(i,:)=En(i,1);
end

figure %% Histogram initial Velocities
histogram(in_vel,40);
hold on
title('Initial Velocities')
xlabel('Velocities');
grid on
hold off

figure %% Histogram final velocities
histogram(fin_vel,40);
hold on
title('Final Velocities')
xlabel('Velocities');
grid on
hold off

figure
plot(time_E,E_Kin,time_E,E_Pot,time_E,E_Tot)
hold on
legend({'E_{Kinetic}','E_{Potential}','E_{Tot}= E_{Kin} + E_{Pot}'},'Location','southwest','FontSize',10)
xlabel('time');
ylabel('Energy');
grid on
hold off