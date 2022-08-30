%% d
I = 1;
Ktot = 1;
Ptot = 1;

kAon = 10;
kAoff = 10;
kIon = 10;
kIoff = 10;
kIcat = 10;
kAcat = 100;

% model equations
dAdt =@(A,AP,I,IK) -kAon*(Ptot - AP)*(A) + kAoff*(AP) + kAcat*(IK); 
dAPdt =@(A,AP,I,IK) kAon*(Ptot - AP)*(A) - kAoff*(AP) - kIcat*(AP);
dIdt = @(A,AP,I,IK) -kIon*(Ktot - IK)*(I) + kIoff*(IK) + kIcat*(AP);
dIKdt = @(A,AP,I,IK) kIon*(Ktot - IK)*(I) - kIoff*(IK) - kAcat*(IK);

dxdt = @(t,x)[dAdt(x(1),x(2),x(3),x(4));
              dAPdt(x(1),x(2),x(3),x(4));
              dIdt(x(1),x(2),x(3),x(4));
              dIKdt(x(1),x(2),x(3),x(4))];

[T, X] = ode45(dxdt, [0,0.6], [0,0,I,0] );

figure; hold on;
plot(T,X(:,1),'-r'); % red for RNA
plot(T,X(:,2),'-', 'color', [0.5 0 1]); % purple for protein
plot(T,X(:,3),'-b')
plot(T,X(:,4),'-g')
ylabel('Molecular concentration (micromolar)')
xlabel('Time (hours)')

%% part e 
nMax = 100;

I = 1;
Ptot = 1;

kAon = 10;
kAoff = 10;
kIon = 10;
kIoff = 10;
kIcat = 10;
kAcat = 100;

Kvalue = [] ;
solution = [];

for Ktot = 1e-3:1e-1:1e2

% model equations
dAdt =@(A,AP,I,IK) -kAon*(Ptot - AP)*(A) + kAoff*(AP) + kAcat*(IK); 
dAPdt =@(A,AP,I,IK) kAon*(Ptot - AP)*(A) - kAoff*(AP) - kIcat*(AP);
dIdt = @(A,AP,I,IK) -kIon*(Ktot - IK)*(I) + kIoff*(IK) + kIcat*(AP);
dIKdt = @(A,AP,I,IK) kIon*(Ktot - IK)*(I) - kIoff*(IK) - kAcat*(IK);

dxdt = @(t,x)[dAdt(x(1),x(2),x(3),x(4));
              dAPdt(x(1),x(2),x(3),x(4));
              dIdt(x(1),x(2),x(3),x(4));
              dIKdt(x(1),x(2),x(3),x(4))];

[T, X] = ode45(dxdt, [0,0.6], [0,0,I,0] );

Kvalue = [Kvalue, Ktot] ;
solution = [solution, X(end,3)];

end
disp(solution)
disp(Kvalue)
figure(2);
hold on
semilogx(Ktot,solution(:,1), '-om')
ylabel('A')
xlabel('Ktot(log)')