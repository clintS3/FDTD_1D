% A 1-D finite-difference time-domain method
% simulation for electromagnetic wave propagation

clear; %test

%Supply some constants
eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(mu0*eps0);
eta = sqrt(mu0/eps0);

% Define S parameter, calculate our range for x, calculate dt.
S = 0.9;
f = 2e9; %Hz
w = 2*pi*f;
lambda0 = c/f; %meters
% number of grid cells per wavelength
N_lambda = 1000;
dx = lambda0/N_lambda; %meters

%make the simulation space
min_x = 0;           %meters
max_x = 3*lambda0;   %meters
x = (min_x:dx:max_x)';
X = x;
mid_x = round(length(x)/2);

% make dt
dt = dx*S/c;

% make the material properties over space
epsilon = eps0.*ones(size(X));
mu = mu0.*ones(size(X));
sigma = zeros(size(X));
sigma_star = zeros(size(X));

% make some material discontinuity
sigma(round(length(x)/2)) = 1;


% calculate C and D parameters from Taflove Ch. 3
Ca = (1-sigma.*dt./2./epsilon)./(1+sigma.*dt./2./epsilon);
Cb = (dt./epsilon./dx)./(1+sigma.*dt./2./epsilon);
Da = (1-sigma_star.*dt./2./mu)./(1+sigma_star.*dt./2./mu);
Db = (dt./mu./dx)./(1+sigma_star.*dt./2./mu);


% prep field component arrays based on the Yee grid cell
Hy = zeros(length(x),2);
Ez = zeros(length(x)-1,2);

% define TFSF locations
TFSF_ind_mx = round(length(x)/4);
TFSF_ind_px = length(x) - round(length(x)/4);

%determine where the wave starts and incident field amplitudes.
x_0 = x(end);
if x_0 == x(1)
    Ez0 = 1;
    Hy0 = -1/377;
    k = w/c;
elseif x_0 == x(end)
    Ez0 = 1;
    Hy0 = 1/377;
    k = -w/c;
end


% make the time-array
min_t = 0;
max_t = 5000*dt;
t = min_t:dt:max_t-dt;
Ez_save = zeros(size(t));

%some parameters for the incident waveform
pulse = 0; %select 1 for pulse
if pulse == 1
    t_0 = 50*dt;%1.25e-9/10;
    p = 15*dt; %pulse_width
    % source_wf = 1*exp(-((t - t_0)/p).^2);
    % figure
    % plot(source_wf)
else
    t_0 = 0;
    p = 0;
end


%prep a figure
fig = figure;
fig.Position  = 1.0e+03*[1.4378    0.0418    1.6344    1.1568];

% run the leap-frog time-update loop
for n = 1:length(t)
 
    %Hy n+1/2 Yee update
    Hy(2:end-1,2) = Da(2:end-1).*Hy(2:end-1,1) + Db(2:end-1).*(Ez(2:end,1) - Ez(1:end-1,1));
    %lower x-TFSF
    Hy(TFSF_ind_mx, 2)   = Hy(TFSF_ind_mx, 2) - dt/mu0/dx*Ez0*source(t(n), t_0, p, x(TFSF_ind_mx), x_0, c, k, w);
    %upper x-TFSF
    Hy(TFSF_ind_px+1, 2) = Hy(TFSF_ind_px+1, 2) + dt/mu0/dx*Ez0*source(t(n), t_0, p, x(TFSF_ind_px), x_0, c, k, w);


    %Ez n+1 Yee update
    Ez(:,2) = Ca(1:end-1).*Ez(:,1) + Cb(1:end-1).*(Hy(2:end,2) - Hy(1:end-1,2));
    %lower x-TFSF
    Ez(TFSF_ind_mx, 2) = Ez(TFSF_ind_mx, 2) - dt/eps0/dx*Hy0*source(t(n)+dt/2, t_0, p, x(TFSF_ind_mx)-dx/2, x_0, c, k, w);
    %upper x-TFSF
    Ez(TFSF_ind_px, 2) = Ez(TFSF_ind_px, 2) + dt/eps0/dx*Hy0*source(t(n)+dt/2, t_0, p, x(TFSF_ind_px)+dx/2, x_0, c, k, w);

    %save E-data at one location
    Ez_save(n) =  Ez(mid_x, 2);
    
    % make a plot in time
    make_plot(fig, x, Ez, Hy, TFSF_ind_mx, TFSF_ind_px, eta, t(n))
    
    %overwrite time-series for next iteration
    [Ez, Hy] = prep_next_iteration(Ez, Hy);


end


%take the fft of the signal and plot it
Ez_save_FFT = fft(Ez_save)./length(Ez_save);
%ample frequency
Fs = 1/dt;
%frequency vector
f_plot = (0:length(t(1:n))-1)./length(t(1:n))*Fs;
f_plot = f_plot(1:length(f_plot)/2);
Ez_save_FFT = Ez_save_FFT(1:length(f_plot));

figure;
plot(f_plot/1e9, 20*log10(abs(Ez_save_FFT)));
xlabel('f (GHz)'); ylabel('dB');
grid on



function source_inc = source(t, t_0, p, x, x_0, c, k, w)
    %be sure to change 'if' pulse statement above!

    %sine waveform
    my_unit_step = ((t - t_0 - k*(x - x_0)/w) > 0);
    source_inc = 1*my_unit_step;%.*sin(w*(t - t_0) - k*(x - x_0));
 
    %pulse waveform
    % source_inc = 1*exp(-((t - t_0 - (x-x_0)/c)/p).^2);
end



function make_plot(fig, x, Ez, Hy, TFSF_ind_mx, TFSF_ind_px, eta, t)
    figure(fig);
    subplot(2,1,1);
    plot(x(1:end-1),  Ez(:,2));
    ylim([-1 1]);
    xlabel('x (meters)');
    title('E_z');
    hold on;
    line([x(TFSF_ind_mx) x(TFSF_ind_mx)], [-1 1]); line([x(TFSF_ind_px) x(TFSF_ind_px)], [-1 1])
    hold off;

    subplot(2,1,2);
    plot(x, eta*Hy(:,2));
    ylim([-1 1]);
    xlabel('x (meters)');
    title('\etaH_y');
    hold on;
    line([x(TFSF_ind_mx) x(TFSF_ind_mx)], [-1 1]); line([x(TFSF_ind_px) x(TFSF_ind_px)], [-1 1])
    hold off;

    sgtitle(sprintf(['Time: ' num2str(round(t/1e-12)) ' picoseconds']));
end


function [Ez, Hy] = prep_next_iteration(Ez, Hy)
    Ez(:,1) = Ez(:,2);
    Ez(:,2) = 0;
    Hy(:,1) = Hy(:,2);
    Hy(:,2) = 0;
end



