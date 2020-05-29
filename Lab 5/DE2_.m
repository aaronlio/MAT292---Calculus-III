function [T, Y] = DE2_(f,t0,tN,y0,y1,h)
N = round(((tN-t0)/h),0);

% Defintion of time vector
T = t0:h:tN;

% Definition of the y(solution) vector
Y = zeros(1,N);
Y(1) = y0;
Y(2) = y0+y1*h;

% Loop to get the solution
for i = 1:N-1
    y_p = (Y(i+1)-Y(i))/h;
    y_dp = f(T(i+2),y_p,Y(i+1));
    Y(i+2) = (h^2)*y_dp+2*Y(i+1)-Y(i);
end
