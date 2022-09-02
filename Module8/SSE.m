function S = SSE(lambdam,theta,alpha)

TimeStep = 0.2;           % Time step (baseline = 0.001)
InitialPopulation = 0.091;    % Initial condition (baseline = 200)
% lambdam = 1;                % Growth rate (baseline = 1)
% theta = 1e3;               % Carrying capacity (baseline = 1000)
% alpha = 2;                  % Exponent (baseline = 2)
endTime = 160;               % End time of simulation (baseline = 10)

% Note: If alpha is 1, this becomes the orginal Logistic Growth model.

% Begin simulation

% Program execution parameters
dT = TimeStep;
tsteps = ceil(endTime/dT) + 1;    % Number of time steps required

% Pre-allocate vectors for faster runtime
Time = zeros(1,tsteps);
 N = zeros(1,tsteps);

N_ex = [0.091, 0.103, 0.136, 0.211, 0.354, 0.415, 0.476, 0.602] ;

% Set initial time and initial population
Time(1) = 0;                % Set initial time to zero
N(1) = InitialPopulation;   % Set initial population

% Euler Method
for n = 2:tsteps
    Time(n) = (n-1)*dT;     

    % Simplified Modified Logistic Growth
        N(n) = N(n-1) + dT*( lambdam*N(n-1)*(1-(N(n-1)/(theta^2))^alpha) ); 

end

S = sum((N(1:100:800) - N_ex).^2) ;

end
