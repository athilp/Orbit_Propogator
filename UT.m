
clear all;

Constants;

UKFConstants;

load('hw02_data.dat');
M = hw02_data;
tVec = M(:,1);
alpha = M(:,2);
delta = M(:,3);

%  tspan = linspace(0, 86400, 2);
%  [Wm,Wc,X,y,Y,Y1,P] = UT(x,Pxx,v, Q,eps, R,params, tspan,S);
% 

P = Pxx;

for k = 1:length(tVec)
        if k == 1
            tspan = 0;
        else
            tspan = linspace(tVec(k-1), tVec(k), 2);
        end
    
        %Time update
    [Wm,Wc,X,y,Y,Y1,Pbar] = UT(x,P,params, tspan,S);
   


    disp('Done particles at t = ')
    tVec(k)

    L = size(X,2);
    m=2;
    zmeas=zeros(m,1);
    Z=zeros(m,L);

    for i = 1:L
        Z(:,i) = H_MeasurementEquation(Y(:,i)) + X(10:11,i);       
        zmeas=zmeas+Wm(i)*Z(:,i);  
    end
    Z1=Z-zmeas(:,ones(1,L));
    Pzz=Z1*diag(Wc)*Z1';  
    Pxz=Y1*diag(Wc)*Z1';  

    z = [alpha(k); delta(k)];
    K_k = Pxz * inv(Pzz);
    x_k = y + K_k*(z-zmeas);
    P_k = Pbar - K_k * Pzz * K_k';

    x = x_k;
    P = P_k;
    disp('Done')

 end
% 
% 

x
P


%Measurement update - at t in tspan, inputs: Wm,Wc,X,y,Y,Y1,P
%[z,P] = Measurmement(Wm,Wc,X,y,Y,Y1,data)








%Transform the particle using the measurement equation


function [Wm,Wc,X,y,Y,Y1,P] = UT(x,Pxx,params, tspan,S)
% n=numel(X0);                                 %numer of states
% 
% alpha=1;                                 %default, tunable
% ki=3-n;                                       %default, tunable
% beta=2;                                     %default, tunable
% lambda=alpha^2*(n+ki)-n;                    %scaling factor
%                                             %scaling factor

                                            %weights
n = S.n;
alpha = S.alpha;
ki = S.ki;
beta = S.beta;
lambda = S.lambda;
Q = S.Q;
R = S.R;
eps = S.eps;
v = S.v;

nx = length(x);
d = length(v);
m = length(eps);

x0 = [x;v;S.eps];
P0 = [Pxx zeros(nx,d+m); zeros(d,nx) Q zeros(d,m); zeros(m, nx+d), R];
X = sigmas(x0,P0,S); 

Wm=[lambda/(n+lambda) 0.5/(n+lambda)+zeros(1,2*n)];  
Wc=[lambda/(n+lambda) + (1-alpha^2+beta) 0.5/(n+lambda)+zeros(1,2*n)]; 

% c=sqrt(n+lambda);
% X=sigmas(X0,P0,c); 
% disp('Done Sigma');%sigma points around x
L=size(X,2);

y=zeros(nx,1);
Y=zeros(nx,L);
% 
% opts = odeset ('RelTol',1e-12, 'AbsTol',1e-30);
% tspan = 0:30:86400;
% 
% STM0 = reshape(eye(6,6), 36,1);
% state0 = [x0;STM0];
% disp('Done');
% [t,state_All] = ode45(@ (t, state) derivs_two_body(state,params), tspan, state0, opts);
% disp('Done');

R = zeros(2,2);
STM0 = eye(6,6);
STM0 = reshape(STM0, 36,1);



opts = odeset ('RelTol',1e-12, 'AbsTol',1e-30);

if tspan == 0
    P = Pxx;
    y = x;
    Y = X(1:6,:);
    Y1= Y-y(:,ones(1,L));
    P=Y1*diag(Wc)*Y1';  
else
    for i = 1:L
    %Initial Conditions
    x0 = X(1,i);
    y0 = X(2,i);
    z0 = X(3,i);
    vx0 = X(4,i);
    vy0 = X(5,i);
    vz0 = X(6,i);

    state0 = [x0;y0;z0;vx0;vy0;vz0];
    params_All = params;
  
    [tVec,state_All] = ode45(@ (t, state) derivs_two_body_state(state,params, [X(7,i); X(8,i); X(9,i)]), tspan, state0, opts);
   
    Y(:,i) = state_All(end,:);       
    y=y+Wm(i)*Y(:,i);       

    end
end

Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1';  
end
