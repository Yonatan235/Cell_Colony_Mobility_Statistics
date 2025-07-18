% === Dynamics of a colony (Active Random Walk) in a fluid ===
clear all

% === Simulation Parameters ===
limt = 10;            % Number of colony sizes to simulate
dt = 0.01;            % Time step
T = 10000.01;         % Total simulation time
n = floor(T/dt);      % Number of time steps
m = 10;               % Number of Monte Carlo simulations
w = 1;                % Noise index (or fixed noise case)

dist = zeros(limt,1);                            % Placeholder for colony displacement
eff_drift_1 = zeros(limt,1);                     % Effective drift
eff_diff_1 = zeros(limt,1);                      % Effective diffusion
Diff_analytical = zeros(limt,m);                 % Analytical diffusion prediction
avg_Diff_analytical = zeros(m,1);                % Average analytical diffusion
D_rr = zeros(limt,1); V_rr = zeros(limt,1);      % Analytical drift and variance
Rot_var_analytical_longtime = zeros(n,limt);     
NewRot_mean_analytical_longtime = zeros(n,limt); 
rot_drift_sd = zeros(n,limt);                    % Std dev of rotational drift
NewRot_var_analytical_longtime = zeros(n,limt);  
NewRot_var_analytical_longtime_sdev = zeros(n,limt);  % Std dev of above
NewRot_var_analytical_longtimesmallnoise = zeros(n,limt);  
Rot_mean_analytical_longtime = zeros(m,limt);    % Mean rotational orientation
Rot_varr_analytical_longtime = zeros(m,limt);    % Rotational variance
serr_drft_1 = zeros(limt,1);                     % Std error of drift

% Colony statistics over time
smx = zeros(n,2,limt);    % Mean position
sfx = zeros(n,2,limt);    % Mean force
svarx = zeros(n,2,limt);  % Position variance
MSD = zeros(n,1,limt);    % Mean squared displacement
serrx = zeros(n,2,limt);  % Std error of position
smth = zeros(n,1,limt);   % Mean angular orientation
svarth = zeros(n,1,limt); % Variance of orientation
serrth = zeros(n,1,limt); % Std error of orientation

% Other placeholders
holderr = zeros(n,1);
speed = zeros(limt,1,m);           % Speed of colony
sumdiff = 0;
sigmat = zeros(limt,1);            % Angular noise SD
ThetaCext = zeros(limt,m,n);       % Store colony angles
XCext = zeros(limt,m,n,2);         % Store colony positions

% Physical properties
a = zeros(limt,1); col_area = zeros(limt,1);
ga_trans_colony = zeros(limt,1);   % Translational drag
ga_rot_colony = zeros(limt,1);     % Rotational drag
sqspeed = zeros(limt,m,n);         % Squared speed
xx = zeros(limt,m,limt);           % Random cell positions (angle-based)
mean_tinitial = zeros(limt,m,limt);% Mean initial cell angular orientation
mean_finitial = zeros(limt,m,limt);% Mean initial cell force

% === Simulation Loop Over Colony Sizes ===
for d = 1:limt
    d
    l = 2*pi*1e-6;                 % Cell + bridge length (m)
    a(d) = l*d/(2*pi);             % Colony radius
    col_area = pi*(a(d))^2;      
    % Initialize arrays per colony size
    XC = zeros(m,n,2);             % Position
    ThetaC = zeros(m,n);           % Orientation
    F_C = zeros(m,n,2);            % Net force
    Torque_C = zeros(m,n);         % Net torque
    t = zeros(n,1);
    holder = zeros(m,n); holder2 = holder;
    holder3 = zeros(2,1);
    F = zeros(d,m,n);              % Force per cell
    thet = zeros(d,m,n);           % Orientation per cell

    % Constants
    K = 1.38064852e-23;            % Boltzmann constant
    temp = 295.15;                 % Temperature (Kelvin)
    ga_trans_colony(d) = (32/3)*a(d)*(1.002e-3);
    ga_rot_colony(d) = (32/3)*(1.002e-3)*(a(d)^3);
    ga_f = 10; ga_t = 10;          % Damping constants
    sigmaf = 3e-13;                % Force noise SD
    sigmat(w) = sqrt(0.002*w);     % Angular noise SD
    vv = l/2;                      % Half cell spacing
    ut = 0.02; uf = 1e-12;       

    for j = 1:m
        if j == 1
            mean_tinitial(1:d,:,d) = ut * randn(d,m);
            mean_finitial(1:d,:,d) = uf + 2*uf * rand(d,m);
            xx(1:d,:,d) = (-vv + 2*vv*rand(d,m))/l;
        end

        % Initialize orientation and force
        for cl = 1:d
            thet(cl,j,1) = mean_tinitial(cl,j,d) + sigmat(w)*randn(1);
            F(cl,j,1) = mean_finitial(cl,j,d) + sigmaf*randn(1);
        end

        ThetaC(j,1) = 0;     % Start angle
        XC(j,1,:) = [0; 0];  % Start at origin

        % === Compute Initial Net Force and Torque ===
        for ko = 1:d
            angle = 2*pi*((ko - 0.5 + xx(ko,j,d))/d) + thet(ko,j,1);
            F_C(j,1,1) = F_C(j,1,1) - F(ko,j,1)*cos(angle);
            F_C(j,1,2) = F_C(j,1,2) - F(ko,j,1)*sin(angle);
            Torque_C(j,1) = Torque_C(j,1) - a(d)*sin(thet(ko,j,1))*F(ko,j,1);
        end

        % === Time Evolution of the System ===
        for i = 1:n-1
            t(i) = (i-1)*dt;
            for k = 1:d
                % OU process for force and orientation
                F(k,j,i+1) = F(k,j,i) - ga_f*(F(k,j,i) - mean_finitial(k,j,d))*dt + sqrt(2*sigmaf^2*ga_f)*randn(1)*sqrt(dt);
                thet(k,j,i+1) = thet(k,j,i) - ga_t*(thet(k,j,i) - mean_tinitial(k,j,d))*dt + sqrt(2*sigmat(w)^2*ga_t)*randn(1)*sqrt(dt);
                Torque_C(j,i) = Torque_C(j,i) - a(d)*sin(thet(k,j,i))*F(k,j,i);
            end

            % Update colony orientation via Euler-Maruyama
            ThetaC(j,i+1) = ThetaC(j,i) + Torque_C(j,i)*dt/ga_rot_colony(d) + sqrt(2*K*temp/ga_rot_colony(d))*randn()*sqrt(dt);

            % Compute net force in the rotated frame
            for k = 1:d
                angle = 2*pi*((k - 0.5 + xx(k,j,d))/d) + thet(k,j,i) + ThetaC(j,i);
                F_C(j,i,1) = F_C(j,i,1) - F(k,j,i)*cos(angle);
                F_C(j,i,2) = F_C(j,i,2) - F(k,j,i)*sin(angle);
            end

            % Euler-Maruyama update for position
            holder3 = XC(j,i,:);
            XC(j,i+1,1) = holder3(1) + F_C(j,i,1)/ga_trans_colony(d)*dt + sqrt(2*K*temp/ga_trans_colony(d))*randn()*sqrt(dt);
            XC(j,i+1,2) = holder3(2) + F_C(j,i,2)/ga_trans_colony(d)*dt + sqrt(2*K*temp/ga_trans_colony(d))*randn()*sqrt(dt);
        end
    end

    % === PostProcessing: Compute Ensemble Statistics ===
    for i = 1:n-1
        smx(i+1,2,w) = mean(XC(:,i+1,2));
        XCext(d,:,i+1,1) = XC(:,i+1,1);
        XCext(d,:,i+1,2) = XC(:,i+1,2);
        svarx(i+1,1,w) = var(XC(:,i+1,1));
        svarx(i+1,2,w) = var(XC(:,i+1,2));
        ThetaCext(w,:,i+1) = ThetaC(:,i+1);

        for j = 1:m
            Fmag = sqrt(F_C(j,i,1)^2 + F_C(j,i,2)^2);
            sqspeed(d,j,i) = (Fmag / ga_trans_colony(d))^2;  % Nondimensional squared speed
        end
    end

    t(n) = n * dt;
end

% === Theoretical Speed and Variance ===
thrt_sqrt_meansqspeed = zeros(limt,1);
thrt_sqrt_stdevsqspeed = zeros(limt,1);
for p = 1:limt
    % Theoretical RMS speed
    thrt_sqrt_meansqspeed(p) = sqrt((p/(ga_trans_colony(p)^2))*((sigmaf^2) + var(mean_finitial(p,:))));
    
    % Standard deviation estimate (analytical)
    thrt_sqrt_stdevsqspeed(p) = sqrt(p^2 * ((2*(uf^2)*bernoulli(2)/(ga_trans_colony(p)^2))^2) + ...
                                      p * (4*(uf^4)*bernoulli(4)/(ga_trans_colony(p)^4)));
end

% for i=1:1
%     errorbar(thrt_sqrt_meansqspeed(2:10),thrt_sqrt_stdevsqspeed(2:10),'b')
%     hold on
% end
