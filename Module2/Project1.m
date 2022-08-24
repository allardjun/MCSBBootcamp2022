

% n - number of days
% x - fraction of caffeinated 

nMax = 21; % max number of days to simulate
c = -0.8;
d = 0.156;

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

% THE MODEL ^
% ------------------------------------------
% THE BEHAVIOR / THE OUTPUT ? 

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
xStart = -2 + (2+2)*rand(1, 100) ;

yStart = -2 + (2+2)*rand(1, 100) ;

figure(4)
plot(xStart, yStart, 'ok');
ylabel('xStart)')
xlabel('yStart')


%% part e

for n=1:nMax
    
    xStart(n+1) =  (xStart(n))^2 - (yStart(n))^2 + c;
    yStart(n+1) = 2*xStart(n)*yStart(n) +  d; 
    
end % finished 

figure(5)
plot(xStart, yStart, '.k');
ylabel('x(n)')
xlabel('y(n)')

figure(5)
hold on
if (xStart(22) >2 & yStart(22)>2 & xStart(22) < -2 & yStart(22) < -2 ) 
    plot(xStart(1), yStart(1), 'or') 
else
    plot(xStart(1),yStart(1),'ob' )
end  

