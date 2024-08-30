x = [-2011.990; -382.065; 6316.376; 5.419783; -5.945319; 1.37398];
Pxx = diag([1; 1; 1 ;10^-6 ;10^-6;10^-6]);


S = struct;

S.v = zeros(3,1);
S.eps = zeros(2,1);


S.n = length(x) + length(S.v) + length(S.eps);                                 %numer of states

S.alpha=1;                                 %default, tunable
S.ki=3-S.n;                                       %default, tunable
S.beta=2;                                     %default, tunable
S.lambda=S.alpha^2*(S.n+S.ki)-S.n;    

S.Q = 10^-16*eye(3,3);
S.R = (10/3600*pi/180)^2*eye(2,2);
%scaling factor
 

