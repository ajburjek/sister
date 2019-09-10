%Two Body Problem
e = 0.01;
a = 6871000; %m, 500 km altitude
i = 90; %deg
w = 0; %deg
RAAN = 0; %deg
f = 0; %deg
M = zeros(1,100); 
range = zeros(1,100); 
M(1) = 0; 
elements = [a, e, i, RAAN, w, M];

%Ground station position
x2 = 228832.435;
y2 =  -5336696.538;
z2 = 3473801.891;

t=1;

mu = 3.9860044e+14; %m^3/s^2

a    = elements(1);%semimajor axis, meters
e    = elements(2);%eccentricity
i    = elements(3);%inclination, degrees
Sini = sind(i); Cosi = cosd(i);
raan = elements(4);%ascending node, degrees
Sinraan = sind(raan); Cosraan = cosd(raan);
w    = elements(5);%aop, degrees
Sinw = sind(w); Cosw = cosd(w);
M    = elements(6);%mean anomaly, degrees

%ekf initials 
%{
pstate = zeros(1,100);
pcov = zeros(1,100);
cstate = zeros(1,100);
ccov = zeros(1,100);
%}
corrected = zeros(1,100);
true = zeros(1,100);

t0 = 0;
n = sqrt(mu/(a^3));
tp = t0 - (M(1)/n);

while t<100
    Sm = sind(M); Cm = cosd(M);
    Q = [((Cosw*Cosraan)-(Sinw*Sinraan*Cosi)) ((-Sinw*Cosraan)-(Cosw*Sinraan*Cosi));...
        ((Cosw*Sinraan)+(Sinw*Cosraan*Cosi)) ((-Sinw*Sinraan)+(Cosw*Cosraan*Cosi));...
        Sinw*Sini Cosw*Sini];

    p = a*(1-e^2);%semilatus rectum
    h = sqrt(mu*p);%mag of angular momentum
    r = p/(1 + (e*cosd(f)));%mag of position vector, meters

    Vr = (h*e/p)*sind(f);%radial velocity, meters/sec
    Vtheta = h/r;%angular velocity, meters/sec

    Xdotstar = (Vr*cosd(f)) - (Vtheta*sind(f));
    Ydotstar = (Vr*sind(f)) + (Vtheta*cosd(f));

    X = Q*[r*cosd(f) Xdotstar; r*sind(f) Ydotstar];

    r = transpose(X(:,1));
    v = transpose(X(:,2));
    X = [r v];

    %Satellite position
    x1 = X(1);
    y1 = X(2);
    z1 = X(3);
    
    %initialState = [0;0]; %FIXME
    %initialState = [r; v]; % get below error w/ this tho i think this should be the initial state
    %FIXME Error using ExtendedKalmanFilter Expected State to be a vector.
    initialState = [x1; y1; z1];

    %Range calculation
    range(t) = sqrt(((x1-x2)^2)+((y1-y2)^2)+((z1-z2)^2));
    %rang = append(rang, range); %FIXME append list
    
    %trackingEKF on Automated Driving System Toolbox
    obj = extendedKalmanFilter(@vdpStateFcn,@vdpMeasurementFcn,initialState);
    obj.ProcessNoise = 1.0;  
    obj.MeasurementNoise = 0.01; % 0.1^2
    
    [PredictedState,PredictedStateCovariance] = predict(obj);
    [CorrectedState,CorrectedStateCovariance] = correct(obj,range(t));

    %{
    pstate(t) = PredictedState;
    pcov(t) = PredictedStateCovariance;
    cstate(t) = CorrectedState;
    ccov(t) = CorrectedStateCovariance;
    %}
    
    corrected(t) = sqrt(((x1-x2)^2)+((y1-y2)^2)+((z1-z2)^2)); %FIXME
    
    true(t) = range(t) - corrected(t);

    %Time step
    
    M(t+1) = n*(t-tp);
    t = t + 1;
end

% plot? 

% i recommend xaxis=time and yaxis=x-pos and repeat for each part of the
% state

% plot lines for (true - see calculation below), corrected state (estimate),
% incoming range measurement (measured)

% true measurement is approximately equal to CorrectedState + the derivative of the measurement equation
% with respect to parameters multiplying the parameter prediction error plus the derivative of the output
% equation with respect to noise multiplying the random part of the noise.
%https://www.coursera.org/lecture/battery-state-of-health/4-5-2-deriving-ekf-method-for-parameter-estimation-Ypuog
% i think for my past project, true = measured-modeled 