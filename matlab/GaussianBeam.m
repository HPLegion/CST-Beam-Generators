%% Script to generate ions for CST simulation
%-----------------
%Generates a beam of particles flying in the direction of the z axis (from neg. to pos.)
%Can simulate different ions and electrons
%Gaussian beams with paraxial beam optics (inaccuracies for big convergence angles)
%Perfectly focusing beams (emittance90 = 0) with exact geometry
%Parallel beam (convangle90 = 0) with adjustable beam size
%This revision still has to be tested and could still containt bugs

%clear all;
warning('not cleared')


%% -------------------------Configuration-------------------------
%general
total_current = 1;  % total current in A
n = 10000;  % number of particles
element_z = 6;  % proton number - 0 for electron -> automatically grabs mass
charge_q = 1;  % put '-1' if using electrons!
accel_voltage = 3.0e4;  % accelerating voltage in V
energy_spread_rel = 0.00;  % Relative enery spread dE_1sig/E

%beam
convangle90_x = rad2deg(0.02/0.0023084507 / 8 /1000);  % x convergence angle for 90% enveloppe in deg
convangle90_y = rad2deg(0.02/0.0023084507 / 8 /1000);  % y convergence angle for 90% enveloppe in deg
emit90_x = 0.02/0.0023084507;  % x emittance in mm mrad (90%)
emit90_y = 0.02/0.0023084507;  % y emittance in mm mrad (90%)
%----
beamsize90_x = 3.;  % only for parallel beam! (convangle90 = 0) in mm
beamsize90_y = 3.;  % only for parallel beam! (convangle90 = 0) in mm

%Positioning
xwaist_z = 264;  % xwaist z coordinate in mm
ywaist_z = 264;  % ywaist z coordinate in mm
startpos_x = 0;  % start position x of the beam in mm
startpos_y = 0;  % start position y of the beam in mm
startpos_z = 0;  % start position z of the beam in mm




%% -------------------------Cast inputs (to SI)-------------------------
%beam
convangle90_x = deg2rad(convangle90_x);  % rad
convangle90_y = deg2rad(convangle90_y);  % rad
emit90_x = emit90_x / 1.e6;  % m rad
emit90_y = emit90_y / 1.e6;  % m rad
%----
beamsize90_x = beamsize90_x / 1000;  % m
beamsize90_y = beamsize90_y / 1000;  % m

%define 1 sigma quantities instead of 90%
convangle_x = convangle90_x / 1.6449;
convangle_y = convangle90_y / 1.6449;
emit_x = emit90_x / 2.7055;
emit_y = emit90_y / 2.7055;
beamsize_x = beamsize90_x / 1.6449;
beamsize_y = beamsize90_y / 1.6449;

%Positioning
xwaist_z = xwaist_z / 1000;  % m
ywaist_z = ywaist_z / 1000;  % m
startpos_x = startpos_x / 1000;  % m
startpos_y = startpos_y / 1000;  % m
startpos_z = startpos_z / 1000;  % m

%% !!!!!!!!Only SI after this point unless explicitly said otherwise!!!!!!!!%%


%% -------------------------Constants-------------------------
C = 299792458;  % speed of light in m/s
Q_E = 1.60217657e-19;  % elementary charge in C
M_E = 9.10938291e-31;  % electron mass in kg
% mean atomic mass for the first one hundred elements in kg
% [0]=electron, [1] = H, [2]=He etc.
ATOMIC_MASSES = 1.660539e-27 * [M_E / 1.660539e-27, 1.008, 4.002602, 6.94, 9.012182, 10.81, 12.011, 14.007, 15.999, 18.9984032, 20.1797,...
    22.98976928, 24.305, 26.9815386, 28.085, 30.973762, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955912, 47.867,...
    50.9415, 51.9961, 54.938045, 55.845, 58.933195, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.9216, 78.96, 79.904,...
    83.798, 85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.96, 98., 101.07, 102.9055, 106.42, 107.8682, 112.411,...
    114.818, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.9054519, 137.327, 138.90547, 140.116, 140.90765, 144.242,...
    145., 150.36, 151.964, 157.25, 158.92535, 162.5, 164.93032, 167.259, 168.93421, 173.054, 174.9668, 178.49,...
    180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592, 204.38, 207.2, 208.9804, 209., 210.,...
    222., 223., 226., 227., 232.03806, 231.03588, 238.02891, 237., 244., 243., 247., 247., 251., 252., 257.];




%% -------------------------Data Generation-------------------------
%Generate data for simulation input file
%Energy, mass, charge, current
mean_energy = abs(charge_q) * accel_voltage *Q_E; % mean energy in J
if energy_spread_rel == 0
    energySI = mean_energy * ones(n,1); % energy w/o spread in J
else
    energySI = normrnd(mean_energy, (mean_energy * energy_spread_rel), n);  % energy in J
end
mSI = ATOMIC_MASSES(element_z+1) * ones(n,1);  % mass in kg --- matlab starts inedxing at 1 thus "+1" to be consistent with python version
mSI = 1.6726219e-27*12* ones(n,1);
warning('mass override for 12-C = 12m_p');
chargeSI = charge_q * Q_E * ones(n,1);  % charge in C
currentSI = total_current / n * ones(n,1);  % current in A

% spatial and angular distribution
% x direction
if convangle_x > 0
    
    % Gaussian mode
    if emit_x > 0
        % calculate covariance matrix at waist
        var_xprime = convangle_x^2; % in rad^2
        var_x = emit_x^2 / var_xprime;  % in m^2
        varmat_x = [var_x, 0.; 0., var_xprime];
        mu_x = [startpos_x, 0.];  % in m and rad
        % evolve back to starting position
        drift_x = [1, (startpos_z - xwaist_z); 0, 1];  % driftmatrix for back evolution
        varmat_x = drift_x * varmat_x * (drift_x');  % transform cov matrix to starting point
        xtemp = mvnrnd(mu_x, varmat_x, n);  % x and x' in m and rad
        xSI = xtemp(:,1);
        xprimeSI = xtemp(:,2);
        
        % perfect focus mode
    elseif emit_x == 0
        xprimeSI = normrnd(0, convangle_x,n,1);  % !1D gauss uses std. dev. NOT variance
        xSI = tan(xprimeSI) * (startpos_z - xwaist_z) + startpos_x;
        
    else
        raise ValueError('illegal emittance value')
    end
    % parallel mode
elseif convangle_x == 0
    xSI = normrnd(startpos_x, beamsize_x, n,1);
    xprimeSI = zeros(n,1);
    
else
    raise ValueError('illegal convergence angle')
end

% y direction
if convangle_y > 0
    
    % Gaussian mode
    if emit_y > 0
        % calculate covariance matrix at waist
        var_yprime = convangle_y^2; % in rad^2
        var_y = emit_y^2 / var_yprime;  % in m^2
        varmat_y = [var_y, 0.; 0., var_yprime];
        mu_y = [startpos_y, 0.];  % in m and rad
        % evolve back to starting position
        drift_y = [1, (startpos_z - ywaist_z); 0,1];  % driftmatrix for back evolution
        varmat_y = drift_y * varmat_y * (drift_y');  % transform cov matrix to starting point
        ytemp = mvnrnd(mu_y, varmat_y, n);  % y and y' in m and rad
        ySI = ytemp(:,1);
        yprimeSI = ytemp(:,2);
        
        % perfect focus mode
    elseif emit_y == 0
        yprimeSI = normrnd(0, convangle_y,n,1);  % !1D gauss uses std. dev. NOT variance
        ySI = tan(yprimeSI) * (startpos_z - ywaist_z) + startpos_y;
        
    else
        raise ValueError('illegal emittance value')
    end
    % parallel mode
elseif convangle_y == 0
    ySI = normrnd(startpos_y, beamsize_y, n,1);
    yprimeSI = zeros(n,1);
    
else
    raise ValueError('illegal convergence angle')
end
% z direction
zSI = startpos_z * ones(n,1);  % z in m

% calculate momenta
pSI = sqrt(energySI.^2 + 2 .* energySI .* mSI * (C^2)) / C;  % relativistic correct momentum in kg*m/s
pzSI = pSI ./ sqrt(1 + tan(xprimeSI).^2 + tan(yprimeSI).^2);  % pz in kg*m/s
pxSI = pzSI .* tan(xprimeSI);  % px in kg*m/s
pySI = pzSI .* tan(yprimeSI);  % py in kg*m/s
pREL = pSI / C ./ mSI;  % momentum in units of beta*gamma
pxREL = pxSI / C ./ mSI;  % px in units of beta*gamma
pyREL = pySI / C ./ mSI;  % py in units of beta*gamma
pzREL = pzSI / C ./ mSI;  % pz in units of beta*gamma


%% -------------------------Saving-------------------------
% Fileheader
header_string = ['%% emit90_x        = %0.5e m rad\n'...
'%% emit90_y        = %0.5e m rad\n'...
'%% convangle90_x   = %0.5e rad\n'...
'%% convangle90_x   = %0.5e rad\n'...
'%% ---\n'...
'%% beamsize90_x    = %0.5e m (only relevant if convangle90_x==0)\n'...
'%% beamsize90_y    = %0.5e m (only relevant if convangle90_y==0)\n'...
'%% ---\n'...
'%% x_start         = %0.5e m\n'...
'%% y_start         = %0.5e m\n'...
'%% z_start         = %0.5e m\n'...
'%% xwaist_z        = %0.5e m\n'...
'%% ywaist_z        = %0.5e m\n'...
'%% N               = %0.5e\n'...
'%% I_total         = %0.5e A\n'...
'%% Z               = %0.5e --- 0 = electron\n'...
'%% Accel. Volt.    = %0.5e V\n'...
'%% --------------------------------------------------\n'...
'%% Uses SI units...\n'...
'%% The momentum (mom) is equivalent to beta * gamma.\n'...
'%% --------------------------------------------------\n'...
'%% Columns: pos_x  pos_y  pos_z  mom_x  mom_y  mom_z  mass  charge  current\n'...
'%% --------------------------------------------------'];
header_string = sprintf(header_string, emit90_x, emit90_y, convangle90_x, convangle90_y, beamsize90_x, beamsize90_y, startpos_x, startpos_y, startpos_z, xwaist_z, ywaist_z, n, total_current, element_z, accel_voltage);
% build export file
timestamp = datestr(datetime('now','TimeZone','local'), 'yyyy-mm-dd-HH-MM-SS');
data = [xSI, ySI, zSI, pxREL, pyREL, pzREL, mSI, chargeSI, currentSI];
% savetxt('ions-output-' + timestamp + '.pid', data, delimiter=' ', fmt='%e', header=header_string, comments='%')
filename = ['ions-output-', timestamp, '.pid'];
% filename = ['12C_gauss',...
%     '_rx_',num2str(beamsize90_x*1000),'_ry_',num2str(beamsize90_y*1000)...
%     '_x_',num2str(startpos_x*1000),'_y_',num2str(startpos_y*1000)...
%     , '.pid'];
warning('filename')
dlmwrite(filename, header_string,'delimiter','');
dlmwrite(filename, data,'-append', 'delimiter',' ','precision','%e')




%% -------------------------Plots-------------------------
rSI = sqrt(xSI.^2 + ySI.^2); % transverse radius for plots
% x phase space
figure
subplot(321)
hist3([xSI, pxREL], [100,100],'EdgeColor','none')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view([0,90]);
title('x phase space')
xlabel('x (m)')
ylabel('px (beta*gamma)')
% y phase space
subplot(322)
hist3([ySI, pyREL], [100,100],'EdgeColor','none')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view([0,90]);
title('y phase space')
xlabel('y (m)')
ylabel('py (beta*gamma)')

% energy distribution
% figure(2)
subplot(323)
hist(energySI/Q_E, 101)
title('Energy Distribution')
xlabel('E (eV)')
subplot(324)
hist3([xSI, ySI], [100,100],'EdgeColor','none')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view([0,90]);
title('real space')
xlabel('x (m)')
ylabel('y (m)')
pbaspect([1 1 1 ])

% momentum
% figure(3)
subplot(325)
plot(rSI, pzREL, '.')
title('longitudinal momentum')
xlabel('r (m)')
ylabel('pz (beta * gamma)')
subplot(326)
plot(rSI, sqrt(pxREL.^2 + pyREL.^2), '.')
title('transverse momentum')
xlabel('r (m)')
ylabel('pt (beta * gamma)')



