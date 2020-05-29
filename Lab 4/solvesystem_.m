function [T, X] = solvesystem_(f,g,t0,tN,x0,h)
N = round(((tN-t0)/h),0);
T = linspace(t0,tN,N);
X = zeros(2,N);
X(1,1) = x0(1,1);
X(2,1) = x0(2,1);

for i = 1:N-1
    x1_SL = f(T(i),X(1,i),X(2,i));
    x1_SR = f(T(i),X(1,i),X(2,i))*h + X(1,i);
    x2_SL = g(T(i),X(1,i),X(2,i));
    x2_SR = g(T(i),X(1,i),X(2,i))*h + X(2,i);
    X(1,i+1) = X(1,i) + (h/2)*(x1_SL + f((T(i)+h),x1_SR,x2_SR));
    X(2,i+1) = X(2,i) + (h/2)*(x2_SL + g((T(i)+h),x1_SR,x2_SR));
end
