% sig = sigma_gr_full(freqs,Ef,gamma,T)
% This function calculates the complex surface conductivity of graphene
% accounting for interband and intraband contributions via Eqns. 3 and 4 
% of "J. Appl. Phys. 103, 64302 (2008)".
% 
% Inputs:
%    freqs - Frequency(ies) at which the conductivity is evaluated in units of Hz
%    Ef    - Graphene Fermi-level in units of eV
%    gamma - Graphene Drude scattering rate in units of Hz
%    T     - Temperature in units of K
% 
% Implementation by Ian Williamson <ian.williamson@utexas.edu>
% based on "G. W. Hanson, Journal of Applied Physics 103, 64302 (2008)"
function sig = sigma_gr_full(freqs,Ef,gamma,T) % C^2/J*sec (S)
e_const = 1.602176565e-19;  % C
h_const = 6.62606957e-34;   % J*sec
h_bar_const = h_const/2/pi;
kB_const = 1.3806488e-23;   % J/K

omega = 2*pi*freqs;         % rad*sec^-1
gamma = 2*pi*gamma;         % rad*sec^-1
Ef_J  = 1.60217657e-19*Ef;
kBT = T*kB_const;

sig_intra = -1j*e_const^2*kBT/pi/h_bar_const^2./(omega-1j*gamma)*(Ef_J/kBT+2*log(exp(-Ef_J/kBT)+1));
sig_inter = -1j*e_const^2/4/pi/h_bar_const*log((2*abs(Ef_J)-(omega-1j*gamma)*h_bar_const)./(2*abs(Ef_J)+(omega-1j*gamma)*h_bar_const));

sig = sig_intra+sig_inter;
end