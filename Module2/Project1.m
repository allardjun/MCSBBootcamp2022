

% n - number of days
% x - fraction of caffeinated

nMax = 21; % max number of days to simulate
c = -0.8005;
d = 0.1560006;

% N = 10; % number of scoops in each jar
x = zeros(1,nMax); % fraction caffeinated
x(1) = 0.1; % initial fraction caffeinated

y = zeros(1,nMax); % fraction caffeinated
y(1) = 0.1; % initial fraction caffeinated

for n=1:nMax

    x(n+1) =  (x(n))^2 - (y(n))^2 + c;
    y(n+1) = 2*x(n)*y(n) +  d;
    %x(n) = x(n-1) - 1/N*x(n-1);

end % finished

figure(1);
plot(x,'ob');
ylabel('x(n)')
xlabel('n')

figure(2)
plot(x, y, 'ok');
ylabel('x(n)')
xlabel('y(n)')

figure(3);
plot(y,'or');
ylabel('y(n)')
xlabel('n')

%% Part d
kMax = 1e5; % change for future parts

xStart = -2 + (2+2)*rand(1, kMax) ;

yStart = -2 + (2+2)*rand(1, kMax) ;

figure(4)
plot(xStart, yStart, 'ok');
ylabel('xStart)')
xlabel('yStart')


%% part e
figure(5);
clf;
hold on

newxMax = 22;


for k= 1:kMax

    x(1) = xStart(k);
    y(1) = yStart(k);

    for j = 1:newxMax
        x(j+1) =  (x(j))^2 - (y(j))^2 + c;
        y(j+1) = 2*x(j)*y(j) +  d;

    end

    if(x(22) < -2 || x(22) > 2 || y(22) < -2 || y(22) > 2 || isnan(x(22)) || isnan(y(22)))

        plot(xStart(k),yStart(k),'.r')
    else
        plot(xStart(k),yStart(k), '.b')
    end
end % finish looping through k

