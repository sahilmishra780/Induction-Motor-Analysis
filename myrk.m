function [tspan,x] = myrk(nameoffunction,tbegin,tfinal,initialconditions,stepsize)

h = (tfinal-tbegin)/(stepsize-1);
tspan = linspace(tbegin,tfinal,stepsize);
x(:,1) = initialconditions;

k_1 = zeros(size(x));
k_2 = zeros(size(x));
k_3 = zeros(size(x));
k_4 = zeros(size(x));

for i=1:(length(tspan)-1)                              
    k_1(:,i) = nameoffunction(tspan(i),x(:,i));
    k_2(:,i) = nameoffunction(tspan(i)+0.5*h,x(:,i)+0.5*h*k_1(:,i));
    k_3(:,i) = nameoffunction((tspan(i)+0.5*h),(x(:,i)+0.5*h*k_2(:,i)));
    k_4(:,i) = nameoffunction((tspan(i)+h),(x(:,i)+k_3(:,i)*h));
    x(:,i+1) = x(:,i) + (1/6)*(k_1(:,i)+2*k_2(:,i)+2*k_3(:,i)+k_4(:,i))*h;  
end

end