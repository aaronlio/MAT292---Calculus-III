%Improved Euler's method
function [T, Y] = IEM_solver(f,t0,tN,y0,h)
N = round(((tN-t0)/h),0);
T = linspace(t0,tN,N);
Y = zeros(1,N);
Y(1) = y0;

for i = 1:N-1
    SL = f(T(i),Y(i));
    SR = f((T(i)+h),(Y(i)+h));
    Y(i+1) = Y(i) + (h/2)*(SL+SR);
end
