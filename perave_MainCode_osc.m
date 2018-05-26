%% PBPL PERiod AVErage 1D FEL simulation code %%
%%% Input deck intended to be compatible with WafFEL, 1D period average, and GENESIS %%%
%% P. Musumeci oscillator version %%
clear all
close all

%% physical constants
global c mu0 e0 me eps0 IA Z0

c = 2.99792458.*10^8;                                            % speed of light
e0 = 1.60217657e-19;                                             % electron charge
me = 9.10938291e-31;                                             % electron mass
eps0 = 8.85418782e-12;                                           % eps_0
mu0 = 1.256637e-6;                                               % mu_0
IA = 17045;                                                      % Alfven current
Z0 = 376.73;                                                     % Impedance of free space         

%% Load the User Determined initial conditions
clear power radfield thetap gammap bunch
param.sigma_t = 125e-15;
recirculate = 0;
tapering_strength = 2;
t1 = tic;
Perave_User_Input_osc;

%% Compute the undulator field
compute_undulator_field_v5h
%% Calculate 1-D FEL parameters
rho1D = 1/param.gamma0*(1/8*param.I/IA*param.K.^2/param.sigmax^2/param.ku^2)^(1/3);
Lgain = param.lambdau/(4*sqrt(3)*pi*rho1D);
Lsat =   param.lambdau/rho1D;
Psat = 1.6*rho1D*param.Ee*param.I;
if param.tapering
    [psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching_factor] = bucket_parameters(param.psir);
    a1 = 2*param.lambda0/param.lambdau*e0*param.E0/me/c^2*sin(param.psir);
    a2 = ((2*param.lambda0/param.lambdau)^1.5)*Z0*e0*param.I*sin(param.psir)^2*capture_fraction*bunching_factor/2/param.A_e/me/c^2;
    pmax_prediction=P0+param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*param.Ee*param.I*capture_fraction;
    etamax = param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*capture_fraction
    bunchlength_rms = param.sigma_t;
    peakcurrent = param.I;
end
%% Run the main integration routine
cavitydetuning = -4;
transmission = 0.6;
sigma_omega = 10;
firstpass =1;
tapering_strength = 2;

for npasses = 1:2
    clear power radfield thetap gammap bunch
    t0 = tic;
    perave_core_v6;
    disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
    perave_postprocessor_v6   
    rad_vs_und(:,npasses) = sum(power,2)*param.lambda0*param.zsep/c;
    rad_vs_beam(:,npasses) = power(end,:);
    Eff(npasses) = Efficiency;
    PL(npasses) = pulselength;
    oldfield(1:param.nslices) =0;
    
    if cavitydetuning>0
    oldfield(1,cavitydetuning+1:cavitydetuning+size(radfield,2)) = radfield(end,:)*transmission;
    else
    oldfield(1,1:1+cavitydetuning+size(radfield,2)) = radfield(end,-cavitydetuning:end)*transmission;    
    end
    pause(0.5)
    %%
    jfreq = 1:param.nslices;
    filter = exp(-(jfreq-param.nslices/2).^2/2/sigma_omega^2);
    figure(8)
    subplot(1,2,1)
    plot(abs(fftshift(fft(oldfield))).*filter);
    subplot(1,2,2)
    filterfield = ifft(ifftshift(fftshift(fft(oldfield) ).*filter));
    plot(abs(filterfield),'g')
    hold on
    plot(abs(oldfield),'r')
    hold off
    pause(0.5)
    oldfield = filterfield;
    firstpass = 0;
end
%% Post-process stuff
figure(100)
plot(max(rad_vs_und),'b')
figure(101)
plot([1:1:param.Nsnap]*param.stepsize,rad_vs_und(:,end),'r')
xlim([0,param.Nsnap*param.stepsize])
figure(102)
plot(PL)


