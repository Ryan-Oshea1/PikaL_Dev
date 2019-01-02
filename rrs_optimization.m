%recover the Rrs  from Trs and Srs measurements at all wavelengths, same as
%the excel spreadsheet, but a bit faster of an optimization process
%function rrs_optimization

%start with the exact same copy of the excel spreadsheet
%create a structure that matches the exact same shape of the RRS_process,
%for ease 
clc
%close allwavelengths_rrs_opt
%%%%% Inputs

% View angle

% Wavelengths
wavelengths_rrs_opt = rrs_sheet(17:115,4); % D 

% Total remote sensing reflectance
trs = rrs_sheet(17:115,5); %E
%trs(1:18) = .015
% Sky Remote Sensing Reflectance
srs = rrs_sheet(17:115,6); %F

%%%% Constants

%aph constants
aph_const = rrs_sheet(17:115,11); %K
aph_coef  = rrs_sheet(17:115,12); %L

%aw
aw  = rrs_sheet(17:115,14); %N
bbw = rrs_sheet(17:115,17); %Q

%gp model constants
g_p_9  = .2;
g_p_10 = .63;
g_p_11 = 2.448;

%gw
gw = .113;

%slope of the CDOM 
sdg = .015;

sub_to_above_1 = .52;

sub_to_above_2 = 1.7;

%%%% Pseudo Constants

%Fresnel Reflectance
fresnel = .022;  %H9
%%
%%%%%
close all
trs_wavelengths = wavelengths_rrs_opt;
trs_wavelengths = 360:5:850; %input camera wavelengths here ZHONGPING, should I get the same values when interpolating to 1 nm?
srs_wavelengths = wavelengths_rrs_opt; %input radiometer wavelengths
% Function for minimization, with 3 variables, P, G, and X with
% minimization function being the error.

[R_s, R_p, fresnel,theta_i_deg,theta_t_deg] = fresnelReflectanceCalculator(1,1.34, 90-40);

%fresnel = .022
%generate initial values for the min function
pgxd = [0,0,0,0];%[0.0136953907906954, 0.00718356450470145, 0.00172280886395222, 0.000587124410745092]; % should not affect IV
[err, pgxd_0] = rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);

%call the function here
minimization_function_rrs_handle = @(pgxd)rrs_optimization_func(pgxd,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);

lower_bounds = [.003,.001,.0001,0]; % ZHONGPING, I set the delta bound to be positive

upper_bounds = [1,1,1,1];%[Inf, Inf, Inf, Inf] ZHONGPING, I set the bounds to be 1 for everything

%silence optimizer output
options = optimoptions('fmincon','Display','notify');
%optimize!
tic
minimization_parameters = fmincon(minimization_function_rrs_handle,pgxd_0,[],[],[],[],lower_bounds,upper_bounds );
toc
%plot the output to check that it was created correctly 
[err, pgxd_0, rrs_mod, rrs_mea, bbp_650] = rrs_optimization_func(minimization_parameters,wavelengths_rrs_opt,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2,trs_wavelengths,srs_wavelengths);

figure; plot(trs_wavelengths,rrs_mea); hold on; plot(trs_wavelengths,rrs_mod); hold on; plot(wavelengths_rrs_opt,trs);hold on; plot(wavelengths_rrs_opt,trs-rrs_mea'); legend('measured','mod','trs', 'difference'); title(['error: ' num2str(err) ' P = ' num2str(minimization_parameters(1)) ' G = ' num2str(minimization_parameters(2)) ' X = ' num2str(minimization_parameters(3)) ' D = ' num2str(minimization_parameters(4)) ' bbp650 = ' num2str(bbp_650)])

%% optimize an entire camera image
%get the real data from the other laptop
%provide the camera wavelengths

%the angle that each spatial pixel is at

%calculate fresnel reflectance for each of these angles

%create a new image out of the old image, solving for each angle in a
%vacuum 

%inside the function, interpret all of the constants to the new camera
%wavelengths, interpret the srs data to camera wavelengths
%do this step first

%could this be made more accurate by assuming spatial uniformity on local
%scales? aka if we oversample can we get a better overall estimate, or does
%averaging do that for us?







%What happens when we apply this correction technique to polarized data?,
%should I calculate Fresnel reflectance differently (assume polarized state
%is totally removed 





