clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We start by supplying our constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(mu0*eps0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We then provide our S parameter, calculate our range for x, and use them
% to calculate delta_t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 0.9;
f = 2e9;%Hz
w = 2*pi*f;
lambda0 = c/f;
N_lambda = 100;
dx = lambda0/N_lambda;%meters

min_x = 0;%meters
max_x = 3*lambda0;%meters


x = (min_x:dx:max_x)';
X = x;

mid_x = round(length(x)/2);

dt = dx*S/c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We provide the permittivity, permeability, and conductivity as a function
% of space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = eps0.*ones(size(X));
mu = mu0.*ones(size(X));
sigma = zeros(size(X));
sigma_star = zeros(size(X));

% sigma(round(length(x)/2)) = 5e7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We calculate the "C" and "D" parameters, as given by Taflove and
% Hagness [2000].  We provide a half-spatial-step version as well
% for simplicity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ca = (1-sigma.*dt./2./eps)./(1+sigma.*dt./2./eps);
Cb = (dt./eps./dx)./(1+sigma.*dt./2./eps);
Da = (1-sigma_star.*dt./2./mu)./(1+sigma_star.*dt./2./mu);
Db = (dt./mu./dx)./(1+sigma_star.*dt./2./mu);

Ca_Ez = Ca(1:length(x)-1);
Cb_Ez = Cb(1:length(x)-1);

Da_Hy = Da;
Db_Hy = Db;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We then set our initial conditions (E=H=J=0 everywhere).  Based
% on the Yee grid, some vectors are longer than others.  We keep 1
% past value of all vectors, except for J.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hy = zeros(length(x),2);
Ez = zeros(length(x)-1,2);

TFSF_ind_mx = 50;
TFSF_ind_px = 200;

x_0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We identify the time steps at which we would like to make plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_t = 0;
max_t = 50/f;
t = min_t:dt:max_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And finally we run the FDTD loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = 0;
output_save_Ez = zeros(round(max_t/dt)+1,5);
eta1 = sqrt(sqrt(-1)*2*pi*f*mu(1)/(sigma(1) + sqrt(-1)*2*pi*f*eps(1)));

fig = figure;
% fig.Position  = [1.8000 41.8000 1920 1080];
fig.Position  = 1.0e+03*[1.4378    0.0418    1.6344    1.1568];

% %Open a video handle
% vidSaveLoc = 'C:\Users\wcs0061\Documents\classes\ELEC3320\2023_Fall_ELEC3320\FDTD_rect_waveguide_2D_TFSF.avi';
% v = VideoWriter(vidSaveLoc);
% open(v);


%k value (Taflove eq. 4.12).
% kx = 2/dx*asin(dx/c/dt*sin(w*dt/2));
pulse = 0;
if pulse == 1
    k = w/c;
    t_0 = 50*dt;%1.25e-9/10;
    p = 15*dt; %pulse_width
    % source_wf = 1*exp(-((t - t_0)/p).^2);
    % figure
    % plot(source_wf)
else
    k = w/c;
    t_0 = 0;
    p = 0;
end






for n = 1:length(t)
    ii = ii + 1;
 
    %Hy Yee update
    Hy(2:end-1,2) = Da_Hy(2:end-1).*Hy(2:end-1,1) + Db_Hy(2:end-1).*(Ez(2:end,1) - Ez(1:end-1,1));
    %lower x-TFSF
    Hy(TFSF_ind_mx, 2)   = Hy(TFSF_ind_mx, 2) - dt/mu0/dx*source(t(n), t_0, p, x(TFSF_ind_mx), x_0, c, k, w);
    %upper x-TFSF
    Hy(TFSF_ind_px+1, 2) = Hy(TFSF_ind_px+1, 2) + dt/mu0/dx*source(t(n), t_0, p, x(TFSF_ind_px), x_0, c, k, w);


    %Ez Yee update
    Ez(:,2) = Ca_Ez.*Ez(:,1) + Cb_Ez.*(Hy(2:end,2) - Hy(1:end-1,2));
    %lower x-TFSF
    Ez(TFSF_ind_mx, 2) = Ez(TFSF_ind_mx, 2) - dt/eps0/dx*(-1/377)*source(t(n)+dt/2, t_0, p, x(TFSF_ind_mx)-dx/2, x_0, c, k, w);
    %upper x-TFSF
    Ez(TFSF_ind_px, 2) = Ez(TFSF_ind_px, 2) + dt/eps0/dx*(-1/377)*source(t(n)+dt/2, t_0, p, x(TFSF_ind_px)+dx/2, x_0, c, k, w);

    Ez_save(ii) =  Ez(mid_x, 2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure-making code below.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(fig);
    subplot(3,1,1);
    plot(x(1:end-1),  Ez(:,2));
    ylim([-1 1]);
    xlabel('x (meters)');
    title('E_z');
    hold on;
    line([x(TFSF_ind_mx) x(TFSF_ind_mx)], [-1 1]); line([x(TFSF_ind_px) x(TFSF_ind_px)], [-1 1])
    hold off;

    subplot(3,1,2);
    plot(x, eta1*Hy(:,2));
    ylim([-1 1]);
    xlabel('x (meters)');
    title('\etaH_y');
    hold on;
    line([x(TFSF_ind_mx) x(TFSF_ind_mx)], [-1 1]); line([x(TFSF_ind_px) x(TFSF_ind_px)], [-1 1])
    hold off;

    sgtitle(sprintf(['Time: ' num2str(round(t(n)/1e-12)) ' picoseconds']));

    drawnow;
    % pause(0.01)  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the arrays/variables for the next iteration.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ez(:,1) = Ez(:,2);
    Ez(:,2) = 0;
    Hy(:,1) = Hy(:,2);
    Hy(:,2) = 0;
end

figure
plot(Ez_save)

Ez_save_FFT = fft(Ez_save)./length(Ez_save);


Fs = 1/dt;

Ff = (0:length(t(1:n))-1)./length(t(1:n))*Fs;


figure
plot(Ff/1e9, abs(Ez_save_FFT))



% close(v);

function source_inc = source(t, t_0, p, x, x_0, c, k, w)

    %be sure to change 'if' pulse statement above!

    my_unit_step = ((t - t_0 - (x - x_0)/c) > 0);
    source_inc = 1*my_unit_step.*sin(w*(t - t_0) - k*(x - x_0));
 
    % source_inc = 1*exp(-((t - t_0 - (x-x_0)/c)/p).^2);


end





