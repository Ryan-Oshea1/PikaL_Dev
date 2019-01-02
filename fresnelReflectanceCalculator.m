%turn this into a funciton where the view angle from nadir is input 
function [R_s, R_p, R_avg,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(n1,n2, view_angle_deg)
%Fresnel Reflectance calculator
n_air = n1;%1.000;

n_seawater = n2;%1.32;

%for normal incidence, the two polarizations are the same (R is the power
%reflection coefficient

R= abs((n_air - n_seawater)/(n_air+n_seawater))^2;


%view_angle_deg = 90 % from nadir
theta_i_deg = 90-view_angle_deg; %in degrees
theta_i_rad = degtorad(theta_i_deg); % convert to radians for use in sine 

theta_t_rad = asin(n_air/n_seawater*sin(theta_i_rad));
theta_t_deg = radtodeg(theta_t_rad);

%for light at an angle
r_s= (n_air * cos(theta_i_rad) - n_seawater*cos(theta_t_rad)) / (n_air * cos(theta_i_rad) + n_seawater * cos(theta_t_rad));
r_p= (n_seawater * cos(theta_i_rad) - n_air*cos(theta_t_rad)) / (n_air * cos(theta_t_rad) + n_seawater * cos(theta_i_rad));

R_s = r_s ^ 2;
R_p = r_p ^ 2;
R_avg = (R_s + R_p)/2;
end