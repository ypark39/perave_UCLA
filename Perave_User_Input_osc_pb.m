% Perave_code_user input
%%%%User entered parameters%%%%%
%% from Duris et al. TESSO paper 
%% Prebuncher Parameters
param.R56buncher=11e-6;
param.phaseshift=1.8;
%% Undulator parameters
param.lambdau = 3.2e-2;                                 % undulator period (m)
param.K = 2.82; %e0*Bfield*/me/c/ku;               % RMS undulator parameter
param.ku = 2.*pi./param.lambdau;                   % undulator wavenumber
lwig=8*param.lambdau;                                                             % Undulator length (m)    
% Tapering options
param.tapering = 0;                                         % tapering (-1 acceleration ; 0 no tapering ; 1 decelation; )    
param.z0 = 0;
param.psir = pi/6;

%% Simulation control options
param.phasespacemovie=0;
param.itdp = 1;
param.prebunching = 0;                                                                  % set to 1 to start from a pre-bunched beam. 

param.changeresphase = 0;
saveoutput=0;
% Set simulation length and # of snapshots
param.delz=1;
param.und_periods = round(lwig/param.lambdau);                         % number of undulator periods to simulate;
param.Nsnap = round(lwig/param.lambdau/param.delz);                % number of snapshots to take over the length of the undulator
param.zsep = 4;                                                              
Nslip=round(param.Nsnap/param.zsep);
param.shotnoise = 1;
param.lambda0 = 266e-9;                                    % seed wavelength (m)
param.k = 2*pi/param.lambda0;                                     % wavenumber in free space
param.nslices =100+4*Nslip+4*round(param.sigma_t/param.zsep/param.lambda0*c);

if(~param.itdp)
    param.nslices = 1;
end
if(param.itdp)
    param.Nsnap = floor(param.und_periods/param.delz);        % Note that the result must be an integer
end
param.stepsize = param.lambdau*param.delz;

%% Electron beam parameters
gamma0 = sqrt(param.k/2/param.ku*(1+param.K^2));param.gamma0=gamma0;          % relativistic gamma factor
param.Np = 256;                                          % # of macroparticles (500-1000 well) 
param.Ee = gamma0*me*c^2/e0;                  % Total e-beam energy (eV)
energyspread = 1*20e-15/param.sigma_t;                                       % Absolute energy spread MeV
param.deltagammarel = energyspread/gamma0/0.511;          % Relative energy spread dgamma/gamma
param.deltagamma = gamma0*param.deltagammarel;
param.bunch = 0.5;                                                   % Initial bunching factor
param.bunchphase = -param.psir-pi/2;                     % Initial bunching phase
param.buncherAmp = 5;


betax=2;
emitx=2e-6;
charge = 300e-12;
%param.sigma_t = 40e-15;
if (param.beamdistribution == 1)
    param.I = charge/sqrt(2*pi)/param.sigma_t   ;           % beam current 
else
    param.I = charge/2/param.sigma_t      ;        % beam current 
end

param.sigmax = sqrt(betax*emitx/gamma0);            % beam radius
param.A_e = 2*pi*param.sigmax^2;                          % beam cross section 
Simulation_temporal_window=param.nslices*param.zsep*param.lambda0/c;

%% radiation parameters
P0 =1e9; param.P0=P0;                                               % Peak input power (W) 
A_mode = param.A_e;                                                     % 1D code. area is same for e_beam and radiation
param.waist = sqrt(A_mode*2/pi);
zr = pi*param.waist^2/param.lambda0;                          % Rayleigh length of seed (m)
param.E0 = sqrt(2*P0/c/eps0/A_mode/2);                        % Assume circular polarization
param.slippage = (param.nslices+112)/2*param.lambda0*param.zsep/c;
param.sigma_l = 2400e-15;

%% Simplifying constants
param.chi2 = e0/me/c^2;
param.chi1=mu0*c/2*param.I/param.A_e;