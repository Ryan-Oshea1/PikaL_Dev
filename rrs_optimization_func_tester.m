function [err,initial_values,rrs_mod,rrs_mea,bbp_650,interpolated_trs,srs,fresnel] =  rrs_optimization_func(pgx,wavelengths,trs,srs,aph_const,aph_coef,aw,bbw,g_p_9,g_p_10,g_p_11,sdg,fresnel,gw,sub_to_above_1,sub_to_above_2, trs_wavelengths, srs_wavelengths)

%calculate everything in the excel spreadsheet
P = pgx(1);
G = pgx(2);
X = pgx(3);
D = pgx(4);

% useable_trs_wavelengths_range = ((trs_wavelengths-400)>0) .* ((trs_wavelengths-1000)<0); %750 to 800 nm (whatever that is)
% utwr = find(useable_trs_wavelengths_range);
% trs_wavelengths = trs_wavelengths(utwr);
% trs = trs(utwr);
%check to see if nan or inf, and interpolate if there is
finite_trs = isfinite(trs);
finite_trs_loc = find(finite_trs);
excel_wavelengths = 400:5:1000;
%size(finite_trs_loc)
interpolated_trs = interp1(trs_wavelengths(finite_trs_loc),trs(finite_trs_loc),excel_wavelengths,'linear','extrap');
trs = interpolated_trs;
trs_wavelengths = excel_wavelengths;


%interpolate everything to the camera wavelengths
srs = interp1(srs_wavelengths,srs,trs_wavelengths);
aph_const = interp1(wavelengths,aph_const,trs_wavelengths);
aph_coef = interp1(wavelengths,aph_coef,trs_wavelengths);
aw = interp1(wavelengths,aw,trs_wavelengths);
bbw = interp1(wavelengths,bbw,trs_wavelengths);


H = trs - srs*fresnel;

H_12_range = ((trs_wavelengths-749)>0) .* ((trs_wavelengths-801)<0); %750 to 800 nm (whatever that is)
H_12_range = find(H_12_range);

H_12 = mean(H(H_12_range)); %750 to 800 nm (whatever that is)

I = H - H_12; %1st guess at Rrs 

aph = (aph_const + aph_coef*log(P))*P;

%find where wavelength 440 is located 
[wavelength_440_min,wavelength_440_loc] = min(abs(trs_wavelengths-440));
[wavelength_555_min,wavelength_555_loc] = min(abs(trs_wavelengths-555));
[wavelength_550_min,wavelength_550_loc] = min(abs(trs_wavelengths-550));

adg = G*exp(sdg*(440-trs_wavelengths));

%this at is not complete ASK ZHONGPING ABOUT THIS
at = aph + aw + adg;

Y = 2*(1-1.2*exp(-.9*I(wavelength_440_loc)/I(wavelength_555_loc))); %440/555

bbp = X*(440./trs_wavelengths).^Y;

rrs_w = gw .* bbw ./(at + bbw +bbp);

bbp_a_bb = bbp./(at + bbw + bbp);

g_p_model = g_p_9 .*(1-g_p_10*exp(-(g_p_11).*bbp_a_bb));

rrs_dp = bbp_a_bb .* g_p_model + rrs_w;

rrs_0_minus = rrs_dp;

sub_to_above = sub_to_above_1./(1-sub_to_above_2.*rrs_0_minus);

rrs_mod = sub_to_above .* rrs_0_minus;

rrs_mea = trs - fresnel .* srs - D;

rrs_mod_rrs_mea = (rrs_mod -rrs_mea).^2;

AA12_range_lower = ((trs_wavelengths-399)>0) .* ((trs_wavelengths-676)<0); %360 to 675 (whatever that is), 430 was start for success
AA12_range_lower = find(AA12_range_lower);

AA12_range_higher = ((trs_wavelengths-749)>0) .* ((trs_wavelengths-851)<0); %750 to 850 nm (whatever that is), 776 was end for success
AA12_range_higher = find(AA12_range_higher);

AA12 = abs(nanmean(rrs_mea(AA12_range_lower) +nanmean(rrs_mea(AA12_range_higher)))); %360 to 675 and  750 to 850

AA13 = (nanmean(rrs_mod_rrs_mea(AA12_range_lower)) + nanmean(rrs_mod_rrs_mea(AA12_range_higher))).^.5;

err = abs(AA13/AA12);

%final plot to compare modeled and measured 
%figure; plot(wavelengths,rrs_mea); hold on; plot(wavelengths,rrs_mod); legend('measured','mod')

initial_values = [.1*(I(wavelength_440_loc)/I(wavelength_555_loc))^(-1.5),.1*(I(wavelength_440_loc)/I(wavelength_555_loc))^(-1.5), 1.5*I(wavelength_550_loc), H_12];

bbp_650 = bbp(59); %ouput just to see if it agrees with the fluorometer bbp values
end