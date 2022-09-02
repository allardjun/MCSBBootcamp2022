%% Part (a) ode45
clear 
close all
%parameters
endtime = 160;
lambda = 0.0163;                % Growth rate (baseline = 1)
theta = 0.847;               % Carrying capacity (baseline = 1000)
alpha = 2.1246;                  % Exponent (baseline = 2)
N = zeros(1,endtime);
N(1) = 0.091;

%equation
f = @(N) lambda.*N*(1 - (N/(theta)^2)^alpha);
dxdt = @(t,x) [f(x(1))];

%Solve 
[X,T] = ode45(dxdt, [0,endtime],N(1)) ;

N_ex = [0.091, 0.103, 0.136, 0.211, 0.354, 0.415, 0.476, 0.602] ;

%plot
figure(1); 
hold on;
plot (X(:,1),T, 'r');
plot([0:20:140],N_ex, 'LineWidth',2)

%% part a euler
clear

% User-defined inputs

TimeStep = 1;           % Time step (baseline = 0.001)
InitialPopulation = 200;    % Initial condition (baseline = 200)
lambdam = 1;                % Growth rate (baseline = 1)
theta = 1e3;               % Carrying capacity (baseline = 1000)
alpha = 2;                  % Exponent (baseline = 2)
endTime = 100;               % End time of simulation (baseline = 10)

% Note: If alpha is 1, this becomes the orginal Logistic Growth model.

% Begin simulation

% Program execution parameters
dT = TimeStep;
tsteps = ceil(endTime/dT) + 1;    % Number of time steps required

% Pre-allocate vectors for faster runtime
Time = zeros(1,tsteps);
N = zeros(1,tsteps);

% Set initial time and initial population
Time(1) = 0;                % Set initial time to zero
N(1) = InitialPopulation;   % Set initial population

% Euler Method
for n = 2:tsteps
    Time(n) = (n-1)*dT;     

    % Simplified Modified Logistic Growth
        N(n) = N(n-1) + dT*( lambdam*N(n-1)*(1-(N(n-1)/theta)^alpha) );    
end

% Plotting
    % Exact Solution to the simplified MLG
       ExactPopulation = theta*((N(1)^alpha)./(N(1)^alpha + (theta^alpha-N(1)^alpha)*exp(-alpha*lambdam*Time))).^(1/alpha);

figure,plot(Time,N,Time,ExactPopulation,'LineWidth',2)
xlabel('Time','FontWeight','bold')
ylabel('Population','FontWeight','bold')
grid on
legend('Numerical solution','Exact solution')
title('Modified Logistic Growth','FontSize',12,'FontWeight','bold')
set(gca,'FontWeight','bold')

%% Part b
clear

% User-defined inputs

TimeStep = 0.2;           % Time step (baseline = 0.001)
InitialPopulation = 0.091;    % Initial condition (baseline = 200)
lambdam = 0.0163;                % Growth rate (baseline = 1)
theta = 0.847;               % Carrying capacity (baseline = 1000)
alpha = 2.1246;                  % Exponent (baseline = 2)
endTime = 160;               % End time of simulation (baseline = 10)

% Note: If alpha is 1, this becomes the orginal Logistic Growth model.

% Begin simulation

% Program execution parameters
dT = TimeStep;
tsteps = ceil(endTime/dT) + 1;    % Number of time steps required

% Pre-allocate vectors for faster runtime
Time = zeros(1,tsteps);
% N = zeros(1,tsteps);

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

figure,plot(Time,N,'LineWidth',2)
hold on
plot([0:20:140],N_ex, 'LineWidth',2)
xlabel('Time','FontWeight','bold')
ylabel('Population','FontWeight','bold')
grid on
legend('Numerical solution','Exact solution', 'Experiment')
title('Modified Logistic Growth','FontSize',12,'FontWeight','bold')
set(gca,'FontWeight','bold')

%Part c
S = sum((N(1:112:end) - N_ex).^2) ;
disp(S);

%part d
% display(SSE(1,1000,2));

%% part e
 x0 = [1e-2, 1e1,2];

 x = fminsearch(@(x)SSE(x(1),x(2),x(3)),x0);

 disp(x);
% 1.2348   -0.3550    2.5605
%1.2079    0.5728    2.5777
% 1.2189   -0.5728    2.6826
% 0.0163    0.8747    2.1246



