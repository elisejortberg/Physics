%Elise Jortberg
%Homework 5 PHYS 4606
%Runge Kutta v Analytic Solutions for 2nd order linear ODE

h = 0.005; %step size
x=[0:h:10];
y = zeros(1, length(x));
g = zeros(1, length(x));
Y= zeros(2, length(x));

for i=1:(length(x)-1)
    k0 = h*[g(i); x(i)*exp(-x(i)) - ((pi^2)/4)*y(i)];
    k1 = h*[(g(i)+k0(2)/2); (x(i)+(h/2))*exp(-x(i)-(h/2)) - ((pi^2)/4)*(y(i)+(k0(1)/2))];
    k2 = h*[(g(i)+k1(2)/2); (x(i)+(h/2))*exp(-x(i)-(h/2)) - ((pi^2)/4)*(y(i)+(k1(1)/2))];
    k3 = h*[(g(i)+k2(2)); (x(i)+(h))*exp(-x(i)-(h)) - ((pi^2)/4)*(y(i)+(k2(1)))];
    
    y(i+1) = y(i)+(1/6)*(k0(1) + 2*k1(1) + 2*k2(1) + k3(1));
    g(i+1) = g(i)+(1/6)*(k0(2) + 2*k1(2) + 2*k2(2) + k3(2));
end
 

%Part II, solve for analytic solution found in #2 and compare
y2= zeros(1, length(x));
denominator = pi*((4+(pi^2))^2);

for i =1:length(x)
    num1 = 4*exp(-x(i));
    num2 = pi*(x(i)*(4+(pi^2))+8) - 2*((pi^2)-4)*exp(x(i))*sin(pi*x(i)/2) - 8*pi*exp(x(i))*cos(pi*x(i)/2);
    numerator = num1*num2;
    y2(i) = numerator/denominator;
end

%Please zoom in on plot to see both. They overlap exactly
figure; 
plot(x, y); hold on,
scatter(x, y2); legend('Runge Kutta', 'Analytic'); title('2nd order linear ODE')