
    N = 100000; %iteration numbers

    %% Define initial variables
    N = 100000; %iteration numbers
    % The modulation of the field
    u1_int = zeros(1,N);
    u2_int = zeros(1,N);
    ui_int = zeros(1,N);
    us_int = zeros(1,N);
    length = 1e-3;% Total length: (unit: m)

    %% Define initial parameters:
    u1_int(1) = sqrt(0.01)/50;
    u2_int(1) = sqrt(0.01)/50;
    us_int(1) =sqrt(0.001)/50;
    ui_int(1) =sqrt(0);

    % wavelength %
    lambda_pump = 1550e-9; % unit: m
    lambda_sig = 1500e-9; % unit: m
    lambda_ider = 1603e-9; % unit: m
    % waveguide losses:
    alpha1 = 3e2; % unit: dB/m
    alpha2 = 3e2; % unit: dB/m
    alphai = 3e2; % unit: dB/m
    alphas = 3e2; % unit: dB/m

    alpha1 = 10^(-alpha1/10); % unit: 1/m
    alpha2 = 10^(-alpha2/10); % unit: 1/m
    alphai = 10^(-alphai/10); % unit: 1/m
    alphas = 10^(-alphas/10); % unit: 1/m

    %% Define / import related parameters:
    % % design parameters:
    % A_pump(1) = ap(1)*exp(phi0_p); % Pump light: (unit: Watt)
    % A_sig(1)  = as(1)*exp(phi0_s); % Signal light:(unit: Watt)
    % A_ider(1) = ai(1)*exp(phi0_i); % Ider light:(unit: Watt)
    % length = 1e-3;% Total length: (unit: m)
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

    delta_beta = abs(beta_i + beta_s - beta_p * 2);

    delta_beta = delta_beta1;
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

    % Calculate the electric field intensity integrated over an infinite area%
    int_sig  = sum(sum((index .* index).*(Ey_sig .* Ey_sig)))*da;
    int_pump = sum(sum((index .* index).*(Ey_pump .* Ey_pump)))*da;
    int_ider = sum(sum((index .* index).*(Ey_ider .* Ey_ider)))*da;

    Gamma1 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_pump*da)))/(int_pump^2);
    Gamma2 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_pump*da)))/(int_pump^2);
    Gammai = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_ider.*conj(temp_Ey_ider).*temp_Ey_ider*da)))/(int_ider^2);
    Gammas = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig ).*temp_Ey_sig .*conj(temp_Ey_sig ).* temp_Ey_sig*da)))/(int_sig ^2);
    gamma1 = 3 * w_p * ng_pump * ng_pump / (4*epsilon0*A0*c^2)*Gamma1;
    gamma2 = 3 * w_p * ng_pump * ng_pump / (4*epsilon0*A0*c^2)*Gamma2;
    gammai = 3 * w_i * ng_ider * ng_ider / (4*epsilon0*A0*c^2)*Gammai;
    gammas = 3 * w_s * ng_sig  * ng_sig  / (4*epsilon0*A0*c^2)*Gammas;

    Gamma12 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_pump*da)))/(int_pump^2);
    Gamma1s = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig ).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_sig *da)))/(int_pump*int_sig );
    Gamma1i = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_ider*da)))/(int_pump*int_ider);
    Gamma21 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_pump*da)))/(int_pump^2);
    Gamma2s = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig ).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_sig *da)))/(int_pump*int_sig );
    Gamma2i = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_pump.*conj(temp_Ey_pump).*temp_Ey_ider*da)))/(int_pump*int_ider);
    Gammai1 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_ider.*conj(temp_Ey_ider).*temp_Ey_pump*da)))/(int_pump*int_ider);
    Gammai2 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_ider.*conj(temp_Ey_ider).*temp_Ey_pump*da)))/(int_pump*int_ider);
    Gammais = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig ).*temp_Ey_ider.*conj(temp_Ey_ider).*temp_Ey_sig *da)))/(int_sig *int_ider);
    Gammas1 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_sig .*conj(temp_Ey_sig ).*temp_Ey_pump*da)))/(int_pump*int_sig );
    Gammas2 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_sig .*conj(temp_Ey_sig ).*temp_Ey_pump*da)))/(int_pump*int_sig );
    Gammasi = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_sig .*conj(temp_Ey_sig ).*temp_Ey_ider*da)))/(int_sig *int_ider);
    gamma12 = 3 * w_p * ng_pump  * ng_pump  / (4*epsilon0*A0*c^2)*Gamma12;
    gamma1s = 3 * w_p * ng_pump  * ng_sig   / (4*epsilon0*A0*c^2)*Gamma1s;
    gamma1i = 3 * w_p * ng_pump  * ng_ider  / (4*epsilon0*A0*c^2)*Gamma1i;
    gamma21 = 3 * w_p * ng_pump  * ng_pump  / (4*epsilon0*A0*c^2)*Gamma21;
    gamma2s = 3 * w_p * ng_pump  * ng_sig   / (4*epsilon0*A0*c^2)*Gamma2s;
    gamma2i = 3 * w_p * ng_pump  * ng_ider  / (4*epsilon0*A0*c^2)*Gamma2i;
    gammai1 = 3 * w_i * ng_ider  * ng_pump  / (4*epsilon0*A0*c^2)*Gammai1;
    gammai2 = 3 * w_i * ng_ider  * ng_pump  / (4*epsilon0*A0*c^2)*Gammai2;
    gammais = 3 * w_i * ng_ider  * ng_sig   / (4*epsilon0*A0*c^2)*Gammais;
    gammas1 = 3 * w_s * ng_sig   * ng_pump  / (4*epsilon0*A0*c^2)*Gammas1;
    gammas2 = 3 * w_s * ng_sig   * ng_pump  / (4*epsilon0*A0*c^2)*Gammas2;
    gammasi = 3 * w_s * ng_sig   * ng_ider  / (4*epsilon0*A0*c^2)*Gammasi;

    Gamma1is2 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_ider.*conj(temp_Ey_sig ).*temp_Ey_pump*da)))/((int_pump*int_pump*int_ider*int_sig)^0.5);
    Gamma2is1 = (A0 * chi3_si * sum(sum(conj(temp_Ey_pump).*temp_Ey_ider.*conj(temp_Ey_sig ).*temp_Ey_pump*da)))/((int_pump*int_pump*int_ider*int_sig)^0.5);
    Gammai12s = (A0 * chi3_si * sum(sum(conj(temp_Ey_ider).*temp_Ey_pump.*conj(temp_Ey_pump ).*temp_Ey_sig*da)))/((int_pump*int_pump*int_ider*int_sig)^0.5);
    Gammas12i = (A0 * chi3_si * sum(sum(conj(temp_Ey_sig ).*temp_Ey_pump.*conj(temp_Ey_pump ).*temp_Ey_ider*da)))/((int_pump*int_pump*int_ider*int_sig)^0.5);
    gamma1is2 = 3 * w_p * sqrt(ng_pump*ng_ider*ng_sig*ng_pump)/(4*epsilon0*A0*c^2)*Gamma1is2;
    gamma2is1 = 3 * w_p * sqrt(ng_pump*ng_ider*ng_sig*ng_pump)/(4*epsilon0*A0*c^2)*Gamma2is1;
    gammai12s = 3 * w_i * sqrt(ng_ider*ng_pump*ng_pump*ng_sig)/(4*epsilon0*A0*c^2)*Gammai12s;
    gammas12i = 3 * w_s * sqrt(ng_sig*ng_pump*ng_pump*ng_ider)/(4*epsilon0*A0*c^2)*Gammas12i;

    for j = 1:N-1

        z = dz*(j-1);
        u1 = u1_int(j);
        u2 = u2_int(j);
        us = us_int(j);
        ui = ui_int(j);

        du1_div_dz = -alpha1/2*u1 + 1i*(gamma1*abs(u1)^2 + 2*(gamma12*abs(u2)^2+gamma1s*abs(us)^2+gamma1i*abs(ui)^2))*u1 + 2*1i*gamma1is2*ui*us*u2'*exp(1i*delta_beta*z);
        du2_div_dz = -alpha2/2*u2 + 1i*(gamma2*abs(u2)^2 + 2*(gamma21*abs(u1)^2+gamma2s*abs(us)^2+gamma2i*abs(ui)^2))*u2 + 2*1i*gamma2is1*ui*us*u1'*exp(1i*delta_beta*z);
        dui_div_dz = -alphai/2*ui + 1i*(gammai*abs(ui)^2 + 2*(gammai1*abs(u1)^2+gammai2*abs(u2)^2+gammais*abs(us)^2))*ui + 2*1i*gammai12s*u1*u2*us'*exp(-1i*delta_beta*z);
        dus_div_dz = -alphas/2*us + 1i*(gammas*abs(us)^2 + 2*(gammas1*abs(u1)^2+gammas2*abs(u2)^2+gammasi*abs(ui)^2))*us + 2*1i*gammas12i*u1*u2*ui'*exp(-1i*delta_beta*z);

        u1_int(j+1) = du1_div_dz * dz +  u1_int(j);
        u2_int(j+1) = du2_div_dz * dz +  u2_int(j);
        ui_int(j+1) = dui_div_dz * dz +  ui_int(j);
        us_int(j+1) = dus_div_dz * dz +  us_int(j);
    end
    
    figure;
    hold on;
    x_axis = (1:N)*dz;
    plot(x_axis,10*log10(abs(u1_int).*abs(u1_int)*1e3));
    plot(x_axis,10*log10(abs(us_int).*abs(us_int)*1e3));
    plot(x_axis,10*log10(abs(ui_int).*abs(ui_int)*1e3));
    legend("P_{pump}","P_{sig}","P_{ilder}")
    xlabel("Propagation distance (m)")
    ylabel("Power (dBm)")
    p_ider = max(abs(ui_int).*abs(ui_int)*1e3);
    p_pump = max(abs(u1_int).*abs(u1_int)*1e3);
    p_sig = max(abs(us_int).*abs(us_int)*1e3);