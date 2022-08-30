%Discrete Logistic growth
clear all
%% part c
nMax = 1000;
r = 0.1;
K = 0.6;

x = zeros(1,nMax);

x(1) = 0.00001;

for n=1:nMax

    x(n+1) = x(n) + r*(1 - (x(n)/K))*x(n);

end % finished loop through n

figure(1);
plot(x,'-ok');
ylabel('rabbit (thousands)')
xlabel('months')

%% part d
r = 2.1 ;
K = 0.6;

x = zeros(1,nMax);

x(1) = 0.1;

for n=1:(nMax +1)

    x(n+1) = x(n) + r*(1 - (x(n)/K))*x(n);

end % finished loop through n

figure(2);
plot(x,'-ob');
ylabel('rabbit (thousands)')
xlabel('months')

%% part e
r = 2.5 ;
K = 0.6;

x = zeros(1,nMax);

x(1) = 0.1;

for n=1:(nMax +1)

    x(n+1) = x(n) + r*(1 - (x(n)/K))*x(n);

end % finished loop through n

figure(3);
plot(x,'-or');
ylabel('rabbit (thousands)')
xlabel('months')

%% part g
figure(4)
clf;

rval = zeros(1, nMax);
for r = 0.1:0.05:2.9
    x = zeros(1, nMax); 
     
    x(1) = 0.1;
    for n=1:nMax

        
        x(n+1) = x(n) + r*(1 - (x(n)/K))*x(n);
        rval(n+1) = r;

    end % finished loop through n
figure(4);

rval = (nMax/2:nMax);
hold on
plot(r, x(nMax/2:nMax),'-om');
ylabel('r')
xlabel('months')

end %finish loop through r



