
%% Clohessy Wiltshire Optimization:
function OrbitRendezvousCurvilinear

%% System Properties
% Propulsion:
S.Tmax = 100; % Maximum Thrust [N]
S.Tmin = 0; % Minimum Thrust [N]
S.ve = 2000; % Emission Velocity [m/s]
S.Re = 6370*1000; % Earth Radius [m]
S.J2 = 1.68262668*10^-3; % J2 Pert

% System:
S.mu = 3.986*10^14; % Gravitational Constant of Earth [m^3/s^2]

% Mass:
S.mdry = 50; % Dry Mass [kg]
S.mwet = 60; % Wet Mass [kg]

% Simulation Time:
T = 10000; %[s]
N = 500; % Number of Steps

% Chaser Orbit Conditions:
% a0 = S.Re + 400*1000; % Semimajor Axis
% e0 = 0; % Eccentricity
% i0 = 110*pi/180; % Inclination
% Om0 = 15*pi/180; % Longitude of the ascending node
% om0 = 0; % Argument of Periapsis
% th0 = -pi/5; % True Anomaly
a0 = S.Re + 400*1000; % Semimajor Axis
e0 = 0.01; % Eccentricity
i0 = 110*pi/180; % Inclination
Om0 = 5*pi/180; % Longitude of the ascending node
om0 = 0; % Argument of Periapsis
th0 = -pi/5; % True Anomaly

% Target Orbit Conditions:
at = S.Re + 840*1000; % Semimajor Axis
et = 0; % Eccentricity
it = 98.6*pi/180; % Inclination
Omt = 0; % Longitude of the ascending node
omt = 0; % Argument of Periapsis
tht = 0; % True Anomaly

S.R = at;
S.n = sqrt(S.mu/S.R^3); % Orbital Rate [1/s]

% State Constraints:
S.rhomin = -600*1000; %[m]
S.rhomax = 600*1000; %[m]

% Scaling:
TU = 1/S.n; % Scaling Time [Seconds per Orbital Periods]
DU = S.R; % Scaling Distance [Meters per Orbital Radius]
MU = S.mwet; % Scaling Mass [kg per Wetmass]

% Apply Scaling:
S.Tmax = S.Tmax*TU^2/(MU*DU);
S.Tmin = S.Tmin*TU^2/(MU*DU);
S.mwet = S.mwet/MU;
S.mdry = S.mdry/MU;
S.ve = S.ve*TU^2/DU;
S.Re = S.Re/DU;
S.R = S.R/DU;
S.n = S.n*TU;
S.mu = S.mu*TU^2/DU^3;
a0 = a0/DU;
at = at/DU;
T = T/TU;
S.rhomin = S.rhomin/DU;
S.rhomax = S.rhomax/DU;

% Convert to Curvilinear Relative Coordinates to get Initial Condition
x0 = [KOE2CURV(S.mu,[a0;e0;i0;Om0;om0;th0],[at;et;it;Omt;omt;tht]); log(S.mwet)];

% Curvilinear System Dynamics(continuous):
S.A = [zeros(3,3) eye(3,3) zeros(3,1);
       3*S.n^2 0 0 0 2*S.n*S.R 0 0;
       0 0 0 -2*S.n/S.R 0 0 0;
       0 0 -S.n^2 0 0 0 0
       0 0 0 0 0 0 0];
   
S.B = [zeros(3,4);
       S.Tmax 0 0 0;
       0 S.Tmax/S.R 0 0;
       0 0 S.Tmax/S.R 0;
       0 0 0 -S.Tmax/S.ve];

%% Optimization:
tf = T;
t = linspace(0,tf,N);
dt = t(2)-t(1);

% Discretization:
% D = expm([S.A S.B; zeros(4,11)]*dt);
% Ad = D(1:7,1:7);
% Bd = D(1:7,8:11);
Ad = expm(S.A*dt);
Bdp = zeros(7,4);
Bdm = zeros(7,4);
for i = linspace(0,dt,20)
    Bdm = Bdm + expm(S.A*(dt-i))*S.B*(dt-i);
    Bdp = Bdp + expm(S.A*(dt-i))*S.B*(i);
end

zmin = log(S.mwet-1/S.ve*S.Tmax*t); %log(max(S.mwet-1/S.ve*S.Tmax*t(k),S.mdry));
zmax = log(S.mwet-1/S.ve*S.Tmin*t);
mu1 = S.Tmin*exp(-zmin);
mu2 = S.Tmax*exp(-zmin);

cvx_begin 

% cvx_solver SeDuMi

variables eta(4,N) x(7,N) s(N)
minimize(sum(s));

subject to
x(:,1) == x0;
x(1:6,N) == zeros(6,1);

for k = 1:N
    if k < N
        x(:,k+1) == Ad*x(:,k) + Bdm*eta(:,k)+Bdp*eta(:,k+1);
%         x(:,k+1) == Ad*x(:,k) + Bd*eta(:,k);
    end
    
    % Mass Constraint:
%     zmin(k) <= x(7,k) <= zmax(k);
    
    % Thrust Constraint
    norm(eta(1:3,k)) <= eta(4,k);
    
    % Upper Thrust Bound
    0 <= eta(4,k) <= mu2(k)*(1-(x(7,k)-zmin(k)));    
    
    % Position Constraints
    S.rhomin <= x(1,k) <= S.rhomax;
    
    % Slack Variable Constraint
    eta(4,k) <= s(k);
    
end
    
cvx_end

%% Plotting:
% Convert Curvilinear Coordinates to Cartesian:
[xt,xc,xcSPHERE,uc] = CURV2CART(S,x,t,[at;et;it;Omt;omt;tht],eta(1:3,:));
xcSPHERE = [xcSPHERE; x(7,:)];
xcSPHERE(2,:) = unwrap(xcSPHERE(2,:));
xcSPHERE(3,:) = unwrap(xcSPHERE(3,:));
uc(4,:) = vecnorm(uc);

% Sequential Convex Program:
% Parameters:
SCP.betash = 2;
SCP.betagr = 2;
SCP.rho0 = 0;
SCP.rho1 = .1;
SCP.rho2 = .7;
SCP.lambda = 100;
SCP.lambda0 = 100;
SCP.lambdaf = 100;
SCP.Nsub = 20;
SCP.eta = .01;
SCP.etau = .5;
SCP.etap = T*.1;
SCP.eta0 = 10^-3;
SCP.eta1 = 10;
SCP.eps = .001;
SCP.imax = 8;

S.pmin = T*.8;
S.pmax = T*1.2;
[xSCP,uSCP,pSCP] = SequentialAlgorithm(S,@Dynamics,SCP,xcSPHERE,uc,T,xcSPHERE(:,1),xcSPHERE(1:6,end),N);

xcSCP = SPHERE2CART(S,xSCP{end});

% Unscale
drho = x(1,:)*DU;
dtheta = x(2,:);
dphi = x(3,:);
drhodot = x(4,:)*DU/TU;
dthetadot = x(5,:)/TU;
dphidot = x(6,:)/TU;
S.R = S.R*DU;
t = t*TU;
T = eta(1:3,:).*exp(x(7,:))*DU*MU/TU^2;
Tmag = vecnorm(T);
m = exp(x(7,:))*MU;
xt(1:3,:) = xt(1:3,:)*DU;
xc(1:3,:) = xc(1:3,:)*DU;
xcSCP = xcSCP*DU;

% figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5)
% plot(t,drho);
% for i = 1:length(xSCP)
%     plot(t,xSCP{i}(1,:) - vecnorm(xt(1:3,:)));
% end
% % plot(t,Y);
% % plot(t,Z);
% xlabel('Time [s]');
% ylabel('Relative Radii [m]');

% figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5)
% plot(t,dtheta*180/pi);
% plot(t,dphi*180/pi);
% xlabel('Time [s]');
% ylabel('Relative Angles [deg]');
% legend('\delta\theta','\delta\phi');
% 
% % figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5)
% % plot(t,Vx);
% % plot(t,Vy);
% % plot(t,Vz);
% % xlabel('Time [s]');
% % ylabel('Velocity [m/s]');
% % legend('X-Velocity','Y-Velocity','Z-Velocity');
% 
% figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5)
% plot(t,Tmag);
% plot(t,T(1,:));
% plot(t,T(2,:));
% plot(t,T(3,:));
% legend('Thrust Magnitude','X-Thrust','Y-Thrust','Z-Thrust');
% xlabel('Time [s]');
% ylabel('Thrust [N]');
% 
% figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5)
% plot(t,m);
% xlabel('Time [s]');
% ylabel('Mass [kg]');

figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5);
plot(t,xcSPHERE(1,:)*DU);
for i = 1:length(xSCP)
    plot(t,xSCP{i}(1,:)*DU);
end
xlabel('Time [s]');
ylabel('$\rho$ [m]');
legend('Initial Guess','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Iteration 5','Iteration 6','Iteration 7','Iteration 8');

figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5);
plot(t,xcSPHERE(2,:));
for i = 1:length(xSCP)
    plot(t,xSCP{i}(2,:));
end
xlabel('Time [s]');
ylabel('$\theta$ [rad]');
legend('Initial Guess','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Iteration 5','Iteration 6','Iteration 7','Iteration 8');

figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5);
plot(t,xcSPHERE(3,:));
for i = 1:length(xSCP)
    plot(t,xSCP{i}(3,:));
end
xlabel('Time [s]');
ylabel('$\phi$ [rad]');
legend('Initial Guess','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Iteration 5','Iteration 6','Iteration 7','Iteration 8');

figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5);
plot(t,m);
for i = 1:length(xSCP)
    plot(t,exp(xSCP{i}(7,:))*MU);
end
xlabel('Time [s]');
ylabel('Mass [kg]');
legend('Initial Guess','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Iteration 5','Iteration 6','Iteration 7','Iteration 8');

figure; hold on; set(gca, 'FontName', 'times','DefaultLineLineWidth', 2,'Box','on','LineWidth',1.5);
plot(t,Tmag);
for i = 1:length(xSCP)
    plot(t,uSCP{i}(4,:)*S.Tmax*DU*MU/TU^2);
end
xlabel('Time [s]');
ylabel('Thrust [N]');
legend('Initial Guess','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Iteration 5','Iteration 6','Iteration 7','Iteration 8');

figure; hold on;

% ChaseAnimate(t,xt(1:3,:),xc(1:3,:),TU);

ChasePlot(xt(1:3,:),xc(1:3,:));
ChasePlot(xt(1:3,:),xcSCP(1:3,:));

end

function xCART = KOE2CART(mu,xKOE)

for j = 1:size(xKOE,2)
    
    a = xKOE(1,j); 
    e = xKOE(2,j); 
    i = xKOE(3,j); 
    Om = xKOE(4,j); 
    om = xKOE(5,j); 
    th = xKOE(6,j);

    p = a*(1-e^2);
    h = sqrt(p*mu);

    r = p/(1+e*cos(th));

    % Determine r, v in pqw frame:
    rpqw = [r*cos(th);r*sin(th);0];
    vpqw = [-mu/h*sin(th);mu/h*(e+cos(th));0];

    % Transfrom pqw to ijk:
    Tpqw2ijk = [cos(Om)*cos(om)-sin(Om)*sin(om)*cos(i), -cos(Om)*sin(om)-sin(Om)*cos(om)*cos(i), sin(Om)*sin(i);
                sin(Om)*cos(om)+cos(Om)*sin(om)*cos(i), -sin(Om)*sin(om)+cos(Om)*cos(om)*cos(i), -cos(Om)*sin(i);
                sin(om)*sin(i), cos(om)*sin(i),cos(i)];
            
    rijk = Tpqw2ijk*rpqw;
    vijk = Tpqw2ijk*vpqw;

    xCART = [rijk;vijk];
    
end

end

function xCART = SPHERE2CART(S,x)

for i = 1:length(x)
    
    xCART(1,i) = x(1,i)*cos(x(2,i))*cos(x(3,i));
    xCART(2,i) = x(1,i)*sin(x(2,i))*cos(x(3,i));
    xCART(3,i) = x(1,i)*sin(x(3,i));
    
end

end

function [xt,xc,xcSPHERE,uc] = CURV2CART(S,x,t,koet,uc)

at = koet(1); 
et = koet(2); 
it = koet(3); 
Omt = koet(4); 
omt = koet(5); 
th0 = koet(6);

thetat = th0 + S.n*t;
pt = at*(1-et^2);
ht = sqrt(pt*S.mu);

for i = 1:length(x)
    
    % To Cartesian:
    rt = pt/(1+et*cos(thetat(i)));
    
    rpqwt = [rt*cos(thetat(i));rt*sin(thetat(i));0];
    vpqwt = [-S.mu/ht*sin(thetat(i));S.mu/ht*(et+cos(thetat(i)));0];
    
    Tpqw2ijkt = [cos(Omt)*cos(omt)-sin(Omt)*sin(omt)*cos(it), -cos(Omt)*sin(omt)-sin(Omt)*cos(omt)*cos(it), sin(Omt)*sin(it);
            sin(Omt)*cos(omt)+cos(Omt)*sin(omt)*cos(it), -sin(Omt)*sin(omt)+cos(Omt)*cos(omt)*cos(it), -cos(Omt)*sin(it);
            sin(omt)*sin(it), cos(omt)*sin(it),cos(it)];
    
    rijkt = Tpqw2ijkt*rpqwt;
    vijkt = Tpqw2ijkt*vpqwt;
        
    xt(:,i) = rijkt;
    
    dx = cos(x(2,i))*cos(x(3,i))*(S.R+x(1,i))-S.R;
    dy = sin(x(2,i))*(S.R+x(1,i));
    dz = sin(x(3,i))*(S.R+x(1,i));
    
    xhat = rijkt/norm(rijkt);
    zhat = cross(rijkt,vijkt)/norm(cross(rijkt,vijkt));
    yhat = cross(zhat,xhat);
    
    xc(1:3,i) = rijkt + xhat*dx + yhat*dy + zhat*dz;
    
    Tv = [cos(x(3,i))*cos(x(2,i)) cos(x(3,i))*sin(x(2,i)) sin(x(3,i)); -sin(x(2,i))/((S.R+x(1,i))*cos(x(3,i))) cos(x(2,i))/((S.R+x(1,i))*cos(x(2,i))) 0;-sin(x(3,i))*cos(x(2,i))/(S.R+x(1,i)) -sin(x(3,i))*sin(x(2,i))/(S.R+x(1,i)) cos(x(3,i))/(S.R+x(1,i))];
    
    uvw = Tv^-1*(x(4:6,i)+[0;1/S.n;0]) - [0;S.R/S.n;0];
    
    xc(4:6,i) = vijkt + xhat*uvw(1) + yhat*uvw(2) + zhat*uvw(3);
    
    % To Sphere:
    rho = norm(xc(1:3,i));
    theta = atan2(xc(2,i),xc(1,i));
    phi = asin(xc(3,i)/rho);
    vrho = dot(xc(1:3,i),xc(4:6,i))/rho;
    vtheta = (xc(1,i)*xc(5,i)-xc(4,i)*xc(2,i))/xc(1,i)^2*rho*cos(phi)*cos(theta)^2;
    vphi = (rho*xc(6,i)-xc(3,i)*vrho)/(rho*cos(phi));
    
    xcSPHERE(:,i) = [rho;theta;phi;vrho;vtheta;vphi];
    
    xchat = xc(1:3,i)/norm(xc(1:3,i));
    zchat = cross(xc(1:3,i),xc(4:6,i))/norm(cross(xc(1:3,i),xc(4:6,i)));
    ychat = cross(zchat,xchat);
    
    uc(1:3,i) = [xchat ychat zchat]'*[xhat yhat zhat]*uc(1:3,i);
    
end

end

function xCURV = KOE2CURV(S,koec,koet)

xc = KOE2CART(S,koec);
xt = KOE2CART(S,koet);

rc = xc(1:3);
vc = xc(4:6);
erc = rc/norm(rc);
ephc = cross(rc,vc)/norm(cross(rc,vc));
ethc = cross(ephc,erc);

rt = xt(1:3);
vt = xt(4:6);
ert = rt/norm(rt);
etht = vt/norm(vt);
epht = cross(rt,vt)/(norm(rt)*norm(vt));

drho = norm(rc) - norm(rt);
dtheta = atan2(dot(erc,etht),dot(erc,ert));
dphi = asin(dot(erc,epht));

drhodot = dot(erc,vc)-dot(ert,vt);
dthetadot = (dot(rc,rt)*(dot(vc,cross(epht,rt))+dot(rc,cross(epht,vt)))-dot(rc,cross(epht,rt))*(dot(vc,rt)+dot(rc,vc)))/(sec(dtheta)*dot(rc,rt))^2;
dphidot = (norm(rc)^2*dot(vc,epht)-dot(rc,epht)*dot(rc,vc))/norm(rc)^3;

xCURV = [drho;dtheta;dphi;drhodot;dthetadot;dphidot];

end


function [xout,uout,pout] = SequentialAlgorithm(S,DynamicsContinuous,SCP,xk0,uk0,p0,xic,xtc,N)

tk = linspace(0,1,N);
xbar = xk0;
ubar = uk0;
pbar = p0;

Nx = size(xbar,1);
Nu = size(ubar,1);

i = 1;

% Convexify Problem, obtain discretized dynamics
[delta,DynamicsDiscrete] = Convexify(S,DynamicsContinuous,SCP,xbar,ubar,pbar,tk);
Jbar = NonLinearCost(xbar,xic,xtc,delta,SCP.lambda,SCP.lambda0,SCP.lambdaf);

f = figure; hold on;
plot(1:length(xk0),xk0(1,:));

while i < SCP.imax
    
    % Solve Problem as a SOCP
    [L,x,u,p,deltaMax] = ConvexSolver(S,DynamicsDiscrete,xbar,ubar,pbar,SCP,xic,xtc,S.pmin,S.pmax,N,Nx,Nu);
    
    plot(1:length(x),x(1,:)); drawnow;
    
    % Check Convergence Criterion
    e = abs((Jbar - L)/abs(Jbar));
    if e <= SCP.eps
        fprintf('SCvx Converged, Max Infeasibility: %d',deltaMax);
        break;
    end
    fprintf('Current Iteration: %d \n Current Error: %d \n Current Max Infeasibility: %d',i,e,deltaMax);
    
    % Convexify Problem, obtain discretized dynamics
    [delta,DynamicsDiscrete] = Convexify(S,DynamicsContinuous,SCP,xbar,ubar,pbar,tk);
    Jstar = NonLinearCost(x,xic,xtc,delta,SCP.lambda,SCP.lambda0,SCP.lambdaf);
    
    % Calculate Convexification Accuracy
    rho = abs((Jbar - Jstar)/(Jbar-L));
    
    % Update trust region
    if rho < SCP.rho0
        % Update eta, dont update trajectory
        SCP.eta = max(SCP.eta0,SCP.eta/SCP.betash);
    elseif rho < SCP.rho1
        % Update eta and trajectory
        SCP.eta = max(SCP.eta,SCP.eta/SCP.betash);
        xbar = x;
        ubar = u;
        pbar = p;
        Jbar = Jstar;
    elseif rho < SCP.rho2
        %Dont update eta, update trajectory
        xbar = x;
        ubar = u;
        pbar = p;
        Jbar = Jstar;
    else
        % Update eta, update trajectory
        SCP.eta = min(SCP.eta1,SCP.betagr*SCP.eta);
        xbar = x;
        ubar = u;
        pbar = p;
        Jbar = Jstar;
    end
    
    xout{i} = xbar;
    uout{i} = ubar;
    pout{i} = pbar;
    i = i + 1;
    
end

close(f);

end

function [L,x,u,p,deltaMax] = ConvexSolver(S,DynamicsDiscrete,xbar,ubar,pbar,SCP,xic,xtc,pmin,pmax,N,Nx,Nu)

Nx0 = size(xic,1);
Nxf = size(xtc,1);

l = [zeros(6,1); -1];
E = eye(Nx,Nx-1);
Eic = eye(Nx0,Nx0-1);
Etc = eye(Nxf,Nxf);

H0 = eye(7,7);
Hf = [1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0];
% Hf = [1 0 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0];

eta = SCP.eta;
etau = SCP.etau;
etap = SCP.etap;
lambda0 = SCP.lambda0;
lambdaf = SCP.lambdaf;
lambda = SCP.lambda;

% Begin Optimization
cvx_begin

variables x(Nx,N) u(Nu,N) p v(Nx-1,N) s(N) vic(Nx0-1) sic vtc(Nxf) stc

minimize(l'*x(1:Nx,N) + lambda*sum(s) + lambda0*sic + lambdaf*stc)

subject to

% Boundary Conditions:
xic == H0*x(:,1) + Eic*vic;
xtc == Hf*x(:,N) + Etc*vtc;

% Slack Variable Constraints:
norm(vic,1) <= sic;
norm(vtc,1) <= stc;

for k = 1:N
    
    % Dynamics:
    if k ~= N
        x(:,k+1) == DynamicsDiscrete.Ad{k}*x(:,k) + DynamicsDiscrete.Bdm{k}*u(:,k) + DynamicsDiscrete.Bdp{k}*u(:,k+1) + DynamicsDiscrete.Sd{k}*p + DynamicsDiscrete.Rd{k} + E*v(:,k);
    end
    
    % Slack Variable Constraints:
    norm(v(:,k),1) <= s(k);
    
    % Parameter Constraints:
    pmin <= p <= pmax;
    
    % Thrust Constraint
    norm(u(1:3,k)) <= u(4,k);
    
    % Upper Thrust Bound
    0 <= u(4,k) <= 1/S.mwet; %exp(-xbar(7,k))*(1-(x(7,k)-xbar(7,k))); 
    
    % State Constraint:
%     S.rhomin <= x(1,:) <= S.rhomax;
    
    % Parameter Trust Region:
    norm(x(:,k)-xbar(:,k),2) <= eta 
    norm(u(:,k)-ubar(:,k),2) <= etau;
    norm(p-pbar,1) <= etap;
end

cvx_end

L = cvx_optval;
deltaMax = max(s);

end

function J = NonLinearCost(xk,xic,xtc,delta,lambda,lambda0,lambdaf)

H0 = eye(7,7);
Hf = [1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0];
J = -xk(7,end) + lambda*sum(delta) + lambda0*norm((H0*xk(:,1)-xic),1) + lambdaf*norm((Hf*xk(:,end)-xtc),1);

end

function [delta,DynamicsMatrices] = Convexify(S,Dynamics,SCP,xk,uk,p,tk)

[Ad,Bdm,Bdp,Sd,Rd,delta] = Propagate(S,@Dynamics,SCP,xk,uk,tk,p);

DynamicsMatrices.Ad = Ad;
DynamicsMatrices.Bdm = Bdm;
DynamicsMatrices.Bdp = Bdp;
DynamicsMatrices.Sd = Sd;
DynamicsMatrices.Rd = Rd;

end

function [Ad,Bdm,Bdp,Sd,Rd,delta] = Propagate(S,Dynamics,SCP,xk,uk,tk,p)

Nx = size(xk,1);
Nu = size(uk,1);
N = size(xk,2);

Feasible = true;

for k = 1:N-1

    phi0 = eye(Nx);
    Bm0 = zeros(Nx,Nu);
    Bp0 = zeros(Nx,Nu);
    S0 = zeros(Nx,1);
    R0 = zeros(Nx,1);
    P0 = Flatten(xk(:,k),phi0,Bm0,Bp0,S0,R0,Nx,Nu);
    
    tspan = linspace(tk(k),tk(k+1),SCP.Nsub);
    P = RK4(@(t,x)Derivatives(t,x,p,S,@Dynamics,tk(k+1),tk(k),uk(:,k+1),uk(:,k),Nx,Nu),tspan,P0);
    
    [xd,Pphi,PBm,PBp,Ps,Pr] = Unflatten(P,Nx,Nu);
    
    delta(k) = norm(xd-xk(:,k+1));
    
%     if delta > SCP.FeasibilityTolerance
%         Feasible = false;
%     end
    
    Ad{k} = Pphi;
    Bdm{k} = Pphi*PBm;
    Bdp{k} = Pphi*PBp;
    Sd{k} = Pphi*Ps;
    Rd{k} = Pphi*Pr;
    

end

end

function [x,Pphi,PBm,PBp,PS,PR] = Unflatten(P,Nx,Nu)

x = P(1:Nx);
Pphi = reshape(P(Nx+1:Nx+Nx^2),[Nx,Nx]);
PBm = reshape(P(Nx+Nx^2+1:Nx+Nx^2+Nx*Nu),[Nx,Nu]);
PBp = reshape(P(Nx+Nx^2+Nx*Nu+1:Nx+Nx^2+2*Nx*Nu),[Nx,Nu]);
PS = P(Nx+Nx^2+2*Nx*Nu+1:2*Nx+Nx^2+2*Nx*Nu);
PR = P(2*Nx+Nx^2+2*Nx*Nu+1:3*Nx+Nx^2+2*Nx*Nu);

end

function P = Flatten(x,Phi,PBm,PBp,PS,PR,Nx,Nu)

P = [x;
     reshape(Phi,[Nx^2,1]);
     reshape(PBm,[Nx*Nu,1]);
     reshape(PBp,[Nx*Nu,1]);
     PS;
     PR];

end

function x = RK4(xdot,tspan,x0)

x = x0;
dt = tspan(2)-tspan(1);

for i = 1:length(tspan)-1
    
    f1 = xdot(tspan(i),x);
    f2 = xdot(tspan(i)+dt/2,x+f1*dt/2);
    f3 = xdot(tspan(i)+dt/2,x+f2*dt/2);
    f4 = xdot(tspan(i)+dt,x+f3*dt);

    x = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

end

function [Pdot] = Derivatives(t,P,p,S,Dynamics,tk1,tk,uk1,uk,Nx,Nu)

[x,Pphi,PBm,PBp,PS,PR] = Unflatten(P,Nx,Nu);

lambdakm = (tk1-t)/(tk1-tk);
lambdakp = 1-lambdakm;

u = lambdakm*uk + lambdakp*uk1;

[f,A,B,S,R] = Dynamics(S,x,u,p);

psi = Pphi^-1;

Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PSdot = psi*S;
PRdot = psi*R;

Pdot = Flatten(Pxdot,Pphidot,PBmdot,PBpdot,PSdot,PRdot,Nx,Nu);

end

function [f,A,B,S,R] = Dynamics(S,x,u,p)

rho = x(1);
theta = x(2);
phi = x(3);
vrho = x(4);
vtheta = x(5);
vphi = x(6);

Tmax = S.Tmax;
ve = S.ve;
mu = S.mu;
J2 = S.J2;
Re = S.Re;

f = p*[vrho;
     vtheta/(rho*cos(phi));
     vphi/rho;
     (vtheta^2+vphi^2)/rho-mu/rho^2-3/4*J2*mu*Re^2/rho^4*(3*cos(2*phi)-1)+Tmax*u(1);
     -vrho*vtheta/rho+vphi*vtheta/rho*tan(phi)+Tmax*u(2);
     -vrho*vphi/rho-vtheta^2/rho*tan(phi)-3/2*J2*mu*Re^2/rho^4*sin(2*phi)+Tmax*u(3);
     -Tmax/ve*u(4)];
     
A =  p*[                                                                                    0, 0,                                                                  0,           1,                              0,                     0, 0
                                                                   -vtheta/(rho^2*cos(phi)), 0,                                 (vtheta*sin(phi))/(rho*cos(phi)^2),           0,               1/(rho*cos(phi)),                     0, 0
                                                                                -vphi/rho^2, 0,                                                                  0,           0,                              0,                 1/rho, 0
      (2*mu)/rho^3 - vtheta^2/rho^2 - vphi^2/rho^2 + (3*J2*Re^2*mu*(3*cos(2*phi) - 1))/rho^5, 0,                                 (9*J2*Re^2*mu*sin(2*phi))/(2*rho^4),           0,                 (2*vtheta)/rho,          (2*vphi)/rho, 0
                                         (vrho*vtheta)/rho^2 - (vphi*vtheta*tan(phi))/rho^2, 0,                                 (vphi*vtheta*(tan(phi)^2 + 1))/rho, -vtheta/rho, (vphi*tan(phi))/rho - vrho/rho, (vtheta*tan(phi))/rho, 0
             (vphi*vrho)/rho^2 + (vtheta^2*tan(phi))/rho^2 + (6*J2*Re^2*mu*sin(2*phi))/rho^5, 0, - (vtheta^2*(tan(phi)^2 + 1))/rho - (3*J2*Re^2*mu*cos(2*phi))/rho^4,   -vphi/rho,       -(2*vtheta*tan(phi))/rho,             -vrho/rho, 0
                                                                                          0, 0,                                                                  0,           0,                              0,                     0, 0];
      
B = p*[   0,    0,    0,        0
        0,    0,    0,        0
        0,    0,    0,        0
     Tmax,    0,    0,        0
        0, Tmax,    0,        0
        0,    0, Tmax,        0
        0,    0,    0, -Tmax/ve];
    

S = f/p;

R = -A*x-B*u;

end

% function xSPHERE = KOE2SPHERE(S,xKOE)
% 
% for j = 1:size(xKOE,2)
%     
%     a = xKOE(1,j); 
%     e = xKOE(2,j); 
%     i = xKOE(3,j); 
%     Om = xKOE(4,j); 
%     om = xKOE(5,j); 
%     th = xKOE(6,j);
% 
%     p = a*(1-e^2);
%     h = sqrt(p*S.mu);
% 
%     r = p/(1+e*cos(th));
% 
%     % Determine r, v in pqw frame:
%     rpqw = [r*cos(th);r*sin(th);0];
%     vpqw = [-S.mu/h*sin(th);S.mu/h*(e+cos(th));0];
% 
%     % Transfrom pqw to ijk:
%     Tpqw2ijk = [cos(Om)*cos(om)-sin(Om)*sin(om)*cos(i), -cos(Om)*sin(om)-sin(Om)*cos(om)*cos(i), sin(Om)*sin(i);
%                 sin(Om)*cos(om)+cos(Om)*sin(om)*cos(i), -sin(Om)*sin(om)+cos(Om)*cos(om)*cos(i), -cos(Om)*sin(i);
%                 sin(om)*sin(i), cos(om)*sin(i),cos(i)];
%             
%     rijk = Tpqw2ijk*rpqw;
%     vijk = Tpqw2ijk*vpqw;
%     
%     rho = norm(rijk,2);
%     theta = atan2(rijk(2),rijk(1));
%     phi = asin(rijk(3)/rho);
%     
%     Tijk2sphere = [cos(theta)*cos(phi) sin(theta)*cos(phi) sin(phi);
%                    -sin(theta) cos(theta) 0;
%                    -cos(theta)*sin(phi) -sin(theta)*sin(phi) cos(phi)];
%     
%     vsphere = Tijk2sphere*vijk;
% 
%     xSPHERE(:,j) = [rho;theta;phi;vsphere];
% 
% end
% 
% end