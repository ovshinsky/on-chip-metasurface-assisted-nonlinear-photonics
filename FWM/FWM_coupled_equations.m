clear all
close all
%% Declare global variables
global vg_pump vg_sig vg_ider k_pump k_ider k_sig

%% Define initial variables
N = 100000; %iteration numbers
A_pump = zeros(1,N); 
A_sig  = zeros(1,N);
A_ider = zeros(1,N);
ap = zeros(1,N); 
as = zeros(1,N); 
ai = zeros(1,N); 

%% Define initial parameters:
ap(1) = 100;
as(1) = 1;
ai(1) = 0;
% initial phase %
phi0_p = 0;
phi0_s = 0;
phi0_i = 0;
% wavelength %
lambda_pump = 1550e-9; % unit: m
lambda_sig = 1500e-9; % unit: m
lambda_ider = 1603e-9; % unit: m
% waveguide losses:
alpha_in_p = 3e2; % unit: dB/m
alpha_in_s = 3e2; % unit: dB/m
alpha_in_i = 3e2; % unit: dB/m
alpha_in_p = 10^(-alpha_in_p/10); % unit: 1/m
alpha_in_s = 10^(-alpha_in_s/10); % unit: 1/m
alpha_in_i = 10^(-alpha_in_i/10); % unit: 1/m

%% Define / import related parameters:
% design parameters:
A_pump(1) = ap(1)*exp(phi0_p); % Pump light: (unit: Watt)
A_sig(1)  = as(1)*exp(phi0_s); % Signal light:(unit: Watt)
A_ider(1) = ai(1)*exp(phi0_i); % Ider light:(unit: Watt)
length = 1e-3;% Total length: (unit: m)
% waveguide geometry:
waveguide_width = 0.5e-6;
waveguide_height = 0.22e-6;
% mode profiles: (FDE span: y:1um z:0.5um; mesh = 100)
load("TE0_1500nm.mat"); % Signal light
Ex_sig = squeeze(mode1_Ex);
Ey_sig = squeeze(mode1_Ey);
Ez_sig = squeeze(mode1_Ez);
neff_sig  = effective_index(1);
ng_sig = 4.553506; % group index
load("TE0_1550nm.mat"); % Pump light
Ex_pump = squeeze(mode1_Ex);
Ey_pump = squeeze(mode1_Ey);
Ez_pump = squeeze(mode1_Ez);
neff_pump = effective_index(1);
ng_pump = 4.762861; % group index
load("TE0_1603nm.mat"); % ider light
Ex_ider = squeeze(mode1_Ex);
Ey_ider = squeeze(mode1_Ey);
Ez_ider = squeeze(mode1_Ez);
neff_ider = effective_index(1);
ng_ider = 4.847830; % group index

% index profile:
index = squeeze(index_x);
mesh_x = squeeze(y); % direction of waveguide width
mesh_y = squeeze(z); % direction of waveguide height

%% Define physical constants:
chi3_si = 3e-18;
c = 299792458;% Light speed: (unit:m/s)
n_si   = max(max(index));
n_sio2 = min(min(index));
n = n_si;
e = 1.602176565e-19; % elementary charge unit:C
N_density = 10e6; % unit: m-3
epsilon0 = 8.854187817e-12;% 真空介电常数 unit: F/M
u0 =  4 * pi *1e-7; % 真空磁导率 unit: V·s/(A·m)

m0 = 9.1093837015e-31; % mass of the electron unit: kg
mce_pie = 0.26*m0;
mch_pie = 0.39*m0;
miu_e = 1417e-4; %Electron mobility unit:m2V-1s-1
miu_h = 470e-4; %Hole Mobility (µh) unit:m2V-1s-1

%% Calculated several parameters:
dz = length/N; % step length
dx = mesh_x(2)-mesh_x(1); % step width
dy = mesh_y(2)-mesh_y(1); % step height
da = dx * dy;
A0 = waveguide_width * waveguide_height;
waveguide_leftpoint = round(0.5*(1e-6-waveguide_width)/dx);
waveguide_rightpoint = round(101-0.5*(1e-6-waveguide_width)/dx);
waveguide_bottompoint = round(0.5*(0.5e-6-waveguide_height)/dy);
waveguide_toppoint =round( 101-0.5*(0.5e-6-waveguide_height)/dy);

vg_sig = c/ng_sig;
vg_pump = c/ng_pump;
vg_ider = c/ng_ider;

w_p = c*2*pi/lambda_pump;
w_s = c*2*pi/lambda_sig ;
w_i = c*2*pi/lambda_ider;

beta_s = neff_sig  * w_s / c;
beta_p = neff_pump * w_p / c;
beta_i = neff_ider * w_i / c;

delta_beta_lin = abs(beta_p - beta_s);
% delta_beta_lin = 0;
% normalization of mode profiles:
normalization_coefficient_pump = sqrt((2*w_p*u0)/beta_s /( (sum(sum(conj(Ex_pump).*Ex_pump+conj(Ey_pump).*Ey_pump+conj(Ez_pump).*Ez_pump))).*da));
normalization_coefficient_sig = sqrt((2*w_s*u0)/beta_p /( (sum(sum(conj(Ex_sig ).*Ex_sig +conj(Ey_sig ).*Ey_sig +conj(Ez_sig ).*Ez_sig ))).*da));
normalization_coefficient_ider = sqrt((2*w_i*u0)/beta_i /( (sum(sum(conj(Ex_ider).*Ex_ider+conj(Ey_ider).*Ey_ider+conj(Ez_ider).*Ez_ider))).*da));

Ex_pump = normalization_coefficient_pump*Ex_pump;
Ey_pump = normalization_coefficient_pump*Ey_pump;
Ez_pump = normalization_coefficient_pump*Ez_pump;
Ex_sig  = normalization_coefficient_sig *Ex_sig ;
Ey_sig  = normalization_coefficient_sig *Ey_sig ;
Ez_sig  = normalization_coefficient_sig *Ez_sig ;
Ex_ider = normalization_coefficient_ider *Ex_ider ;
Ey_ider = normalization_coefficient_ider *Ey_ider ;
Ez_ider = normalization_coefficient_ider *Ez_ider ;

% Confinement factor k:
temp_Ey_pump = Ey_pump(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ex_pump = Ex_pump(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ez_pump = Ez_pump(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);

temp_Ey_sig  = Ey_sig (waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ex_sig  = Ex_sig(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ez_sig  = Ez_sig(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);

temp_Ey_ider = Ey_ider (waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ex_ider = Ex_ider(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
temp_Ez_ider = Ez_ider(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);

k_sig = (n^2 * sum(sum(temp_Ey_sig .*temp_Ey_sig + temp_Ex_sig .*temp_Ex_sig + temp_Ez_sig .*temp_Ez_sig)))  /  (sum(sum((index .* index).*(Ey_sig .* Ey_sig + Ex_sig .* Ex_sig +Ez_sig .* Ez_sig ))));
k_pump = (n^2 * sum(sum(temp_Ey_pump .*temp_Ey_pump + temp_Ex_pump .*temp_Ex_pump + temp_Ez_pump .*temp_Ez_pump)))  /  (sum(sum((index .* index).*(Ey_pump .* Ey_pump + Ex_pump .* Ex_pump +Ez_pump .* Ez_pump ))));
% k_sig  = (n^2 * sum(sum(temp_Ey_sig  .*temp_Ey_sig )))  /  (sum(sum((index .* index).*(Ey_sig  .* Ey_sig))));
k_ider = (n^2 * sum(sum(temp_Ey_ider .*temp_Ey_ider + temp_Ex_ider .*temp_Ex_ider + temp_Ez_ider .*temp_Ez_ider)))  /  (sum(sum((index .* index).*(Ey_ider .* Ey_ider + Ex_ider .* Ex_ider +Ez_ider .* Ez_ider ))));


% calculate absorptive loss via FCA alpha_FC%
alpha_FC_p = e^3 * N_density / (epsilon0 * c * n * w_p ^2 ) *(1/(miu_e * mce_pie) + 1 / miu_h*(mch_pie)^2);
alpha_FC_s = e^3 * N_density / (epsilon0 * c * n * w_s ^2 ) *(1/(miu_e * mce_pie) + 1 / miu_h*(mch_pie)^2);
alpha_FC_i = e^3 * N_density / (epsilon0 * c * n * w_i ^2 ) *(1/(miu_e * mce_pie) + 1 / miu_h*(mch_pie)^2);

%calculate the free-carrier-induced change in the refractive index%
delta_nFC_p = -e^2/(2*epsilon0*n*w_p^2)*(N_density/mce_pie + N_density^0.8/mch_pie);
delta_nFC_s = -e^2/(2*epsilon0*n*w_s^2)*(N_density/mce_pie + N_density^0.8/mch_pie);
delta_nFC_i = -e^2/(2*epsilon0*n*w_i^2)*(N_density/mce_pie + N_density^0.8/mch_pie);
% 
% delta_nFC_p = 0;
% delta_nFC_s = 0;
% delta_nFC_i = 0;
    
% Calculate the electric field intensity integrated over an infinite area%
int_sig  = sum(sum((index .* index).*(Ey_sig .* Ey_sig)))*da;
int_pump = sum(sum((index .* index).*(Ey_pump .* Ey_pump)))*da;
int_ider = sum(sum((index .* index).*(Ey_ider .* Ey_ider)))*da;
% Calculate the effective susceptibility%
Gamma_pppp = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_pump.*conj(temp_Ey_pump).*(temp_Ey_pump*da))))/(int_pump * int_pump * int_pump * int_pump);
Gamma_spps = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig).*temp_Ey_pump.*conj(temp_Ey_pump).*(temp_Ey_sig*da))))/(int_sig * int_pump * int_pump * int_sig);
Gamma_spip = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig).*temp_Ey_pump.*conj(temp_Ey_ider).*(temp_Ey_pump*da))))/(int_sig * int_pump * int_ider * int_pump);
Gamma_ippi = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_pump.*conj(temp_Ey_pump).*(temp_Ey_ider*da))))/(int_ider * int_pump * int_pump * int_ider);
Gamma_ipsp = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_pump.*conj(temp_Ey_sig).*(temp_Ey_pump*da))))/(int_ider * int_pump * int_pump * int_ider);
%% Solve coupled equations
for j = 1:N-1
    % Current position:
    z = dz*(j-1);
    
    % Amplitude coupled relations:
    dap_div_dz = -c * kp(z) / (2 * n * vg_p(z)) * (alpha_in_p + alpha_FC_p) * ap(j) + ...
                1i * w_p * kp(z) / (n * vg_p(z)) * delta_nFC_p(1) * ap(j) + ...
                1i * 3 * w_p / (4 * epsilon0 * A0(1) * vg_p(z)) * (Gamma_pppp(1) / vg_p(z) * abs(ap(j))^2) * ap(j);

    das_div_dz = -c * ks(z) / (2 * n * vg_s(z)) * (alpha_in_s + alpha_FC_s(1)) * as(j) + ...
                1i * w_s * ks(z) /(n * vg_s(z)) * delta_nFC_s(1) * as(j) + ...
                1i * 3 * w_s / (4 * epsilon0 * A0(1) * vg_s(z)) * (2 * Gamma_spps(1) / vg_p(z) * abs(ap(j))^2) * as(j) + ...  
                1i * 3 * w_i / (4 * epsilon0 * A0(1) * vg_i(z)) * Gamma_spip(1) / vg_p(z) * ap(j)^2 * ai(j)' * exp(1i * delta_beta_lin(1));

    dai_div_dz = -c * ki(z) /(2 * n * vg_i(z)) * (alpha_in_i + alpha_FC_i(1)) * ai(j) + ...
                1i * w_i * ki(z) / (n * vg_i(z)) * delta_nFC_i(1) * ai(j) + ...
                1i * 3 * w_i / (4 * epsilon0 * A0(1) * vg_i(z)) * (2 * Gamma_ippi(1) /vg_p(z) * abs(ap(j))^2) * ai(j) + ...
                1i * 3 * w_s / (4 * epsilon0 * A0(1) * vg_s(z)) * Gamma_ipsp(1) / vg_p(z) * ap(j)^2 * as(j)' * exp(1i * delta_beta_lin(1));

    % Three coupled equations:
    A_pump(j+1) =  A_pump(j) + dap_div_dz * exp(1i * beta_p * z) * dz + 1i * beta_p * A_pump(j) * dz;
    A_sig(j+1)  =  A_sig(j) + das_div_dz * exp(1i * beta_s * z) * dz + 1i * beta_s * A_sig(j)  * dz;
    A_ider(j+1) =  A_ider(j) + dai_div_dz * exp(1i * beta_i * z) * dz + 1i * beta_i * A_ider(j) * dz;

    % Change of ap,ai,as
    ap(j+1) = ap(j)+dap_div_dz * dz;
    as(j+1) = as(j)+dap_div_dz * dz;
    ai(j+1) = ai(j)+dap_div_dz * dz;

    
end

%% Plot:
figure;
hold on;
x_axis = (1:N)*dz;
plot(x_axis,abs(A_pump));
plot(x_axis,abs(A_sig ));
plot(x_axis,abs(A_ider));
legend("A_{pump}","A_{sig}","A_{ider}")
xlabel("propagation distance(m)")
ylabel("amplitude(v/m)")
figure;
hold on;
x_axis = (1:N)*dz;
plot(x_axis,abs(ap));
plot(x_axis,abs(as ));
plot(x_axis,abs(ai));
legend("A_{pump}","A_{sig}","A_{ider}")
xlabel("propagation distance(m)")
ylabel("amplitude(v/m)")

figure;
hold on;
x_axis = (1:N)*dz;
plot(x_axis,log(abs(ap)));
plot(x_axis,log(abs(as )));
plot(x_axis,log(abs(ai)));
legend("A_{pump}","A_{sig}","A_{ider}")
xlabel("propagation distance(m)")
ylabel("log(amplitude(v/m))")


%% Calculated several used parameters (function):
% % Confinement factor k:
    % temp_Ey_pump = Ey_pump(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
    % temp_Ey_sig  = Ey_sig (waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
    % temp_Ey_ider = Ey_ider(waveguide_leftpoint:waveguide_rightpoint, waveguide_bottompoint:waveguide_toppoint);
    % 
    % k_pump = (n^2 * sum(sum(temp_Ey_pump .*temp_Ey_pump)))  /  (sum(sum((index .* index).*(Ey_pump .* Ey_pump))))
    % k_sig  = (n^2 * sum(sum(temp_Ey_sig  .*temp_Ey_sig )))  /  (sum(sum((index .* index).*(Ey_sig  .* Ey_pump))))
    % k_ider = (n^2 * sum(sum(temp_Ey_ider .*temp_Ey_ider)))  /  (sum(sum((index .* index).*(Ey_ider .* Ey_ider))))

function temp = kp(~)
global k_pump
temp = k_pump;
end 

function temp = ks(~)
global k_sig
temp = k_sig;
end 
function temp = ki(~)
global k_ider
temp = k_ider;
end 
function temp = vg_p(~)
global vg_pump
temp = vg_pump;
end 
function temp = vg_i(~)
global vg_ider
temp = vg_ider;
end 
function temp = vg_s(~)
global vg_sig
temp = vg_sig;
end 