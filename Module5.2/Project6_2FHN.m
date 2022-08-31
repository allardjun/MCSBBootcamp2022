% model parameters
eps = 0.08;
a = 1;
b = 0.2;

%I
I0 = 1.0;
tStart = 40;
tStop = 47;
I = @(t) I0*(t>tStart).*(t<tStop);

% model definition
f = @(v,w) v - 1/3*v.^3 - w;
g = @(v,w) eps*(v + a -b*w);

h = @(v,w,t) f(v,w) + I(t); %new f with I included 


%% single cell
dxdt =@ (t,x) [h(x(1),x(2),t); g(x(1),x(2));];

% solve!
[T,X] = ode45(dxdt,[0,100], [-1.5,-0.5]);



% figure(4); hold on;
% set(gca, 'xlim', [-2.5, 2.5], 'ylim', [-2.5,2.5])
% ylabel('w');
% xlabel('v')

uArray = linspace(-2.5, 2.5,32);
wArray = linspace(-2.5, 2.52,32);

% [uMesh,wMesh] = meshgrid(uArray, wArray);
% 
% % the Matlab plot command for a field of arrows is:
% quiver(uMesh, wMesh, f(uMesh, wMesh), g(uMesh,wMesh), 0.5)

figure(5); hold on;
plot (T,X(:,1), 'r');
plot (T,X(:,2),'-', 'color',[0.4940 0.1840 0.5560]);

%% part 4 look in to use of mod().
% nMax = 10;
% Varray = zeros(1,nMax);
% Warray = zeros(1,nMax);
% 
% for i = 1:10
%     % model parameters
% eps = 0.08;
% a = 1;
% b = 0.2;
% D = 0.9;
% 
% %I
% I0 = 1.0;
% tStart = 40;
% tStop = 47;
% I = @(t) I0*(t>tStart).*(t<tStop);
% 
% % model definition
% f = @(v,w) v(i) - 1/3*v(i).^3 - w(i) + D*(v(i-1)-2*v(i)+v(i+1));
% g = @(v,w) eps*(v + a -b*w);
% 
% h = @(v,w,t) f(v,w) + I(t); %new f with I included 
% end
% 
% %single cell
% dxdt =@ (t,x) [h(x(1),x(2),t); g(x(1),x(2));];
% %solve
% [T,X] = ode45(dxdt,[0,100],rand(20,1));

% model parameters
eps = 0.08;
a = 1;
b = 0.2;
D = 0.9;

%I
I0 = 1.0;
tStart = 40;
tStop = 47;
I = @(t) I0*(t>tStart).*(t<tStop);

%equations
% model definition
f(1) = @(v,w) v(1) - 1/3*v(1).^3 - w(1) ;
g(1) = @(v,w) eps*(v + a -b*w);

fnew = @(v,w,t) f(v,w) + I(t) - D*(v(10) - 2*v(1) +v(2)); %new f with I included 





%% single cell
dxdt =@ (t,x) [fnew(x(1),x(2),t); g(x(1),x(2));];

% solve!
[T,X] = ode45(dxdt,[0,100], [-1.5,-0.5]);




