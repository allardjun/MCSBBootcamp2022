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
f = @(v,w) v - 1/3*v.^3 - w ;
g = @(v,w) eps*(v + a -b*w);

fnew{1} = @(v,w,t) f(v,w) + I(t) - D*(v(10) - 2*v(1) +v(2)); %new f with I included 
fnew{2} = @(v,w,t) f(v,w) + I(t) - D*(v(1) - 2*v(2) +v(3)); %new f with I included
fnew{3} = @(v,w,t) f(v,w) + I(t) - D*(v(2) - 2*v(3) +v(4)); %new f with I included
fnew{4} = @(v,w,t) f(v,w) + I(t) - D*(v(3) - 2*v(4) +v(5)); %new f with I included
fnew{5} = @(v,w,t) f(v,w) + I(t) - D*(v(4) - 2*v(5) +v(6)); %new f with I included
fnew{6} = @(v,w,t) f(v,w) + I(t) - D*(v(5) - 2*v(6) +v(7)); %new f with I included
fnew{7} = @(v,w,t) f(v,w) + I(t) - D*(v(6) - 2*v(7) +v(8)); %new f with I included
fnew{8} = @(v,w,t) f(v,w) + I(t) - D*(v(7) - 2*v(8) +v(9)); %new f with I included
fnew{9} = @(v,w,t) f(v,w) + I(t) - D*(v(8) - 2*v(9) +v(10)); %new f with I included
fnew{10} = @(v,w,t) f(v,w) + I(t) - D*(v(9) - 2*v(10) +v(1)); %new f with I included

 g{1} = g;
 g{2} = g;
 g{3} = g;
 g{4} = g;
 g{5} = g;
g{6} = g;
g{7} = g;
g{8} = g;
g{9} = g;
g{10} = g;

F = [fnew(1); fnew(2); fnew(3); fnew(4); fnew(5); fnew(6); fnew(7); fnew(8); fnew(9); fnew(10)];
% G = [g(1); g(2); g(3); g(4); g(5); g(6); g(7); g(8); g(9); g(10)];


%% single cell
dxdt =@ (t,x) [rand(20,1)];

% solve!
[T,X] = ode45(dxdt,[0,100], rand(20,1));




