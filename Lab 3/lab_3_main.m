%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
clear; close all; clc;
f_1 = @(t,y) y*(tan(t))+(sin(t));
t0 = 0;
tN = pi;
y0 = -1/2;
soln_1 = ode45(f_1,[t0, tN], y0);
[IEM_X, IEM_Y] = IEM_solver(f_1,t0,tN,y0,0.01);
subplot(2,2,1);
plot(soln_1.x, soln_1.y, IEM_X, IEM_Y);
legend('ODE45', 'Improved Euler', 'Location', 'Best');
title("y' = y*(tan(t))+(sin(t))");
ylabel('y');
xlabel('t');
% Major differences:
% At pi/2, there is small jump in the Improved Euler's Method (IEM). When
% appraoching from the left, IEM gives a large dy/dx. As a result, the IEM
% solution spikes higher than the ODE45 solution, which also creates a
% large deviation on the right side of pi/2.
% This problem is not seen in the ODE45 solution since it adapts its step
% sizes.

% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
f_1 = @(t,y) 1/(y^2);
t0 = 1;
tN = 10;
y0 = 1;
soln_1 = ode45(f_1,[t0, tN], y0);
[IEM_X, IEM_Y] = IEM_solver(f_1,t0,tN,y0,0.01);
subplot(2,2,2);
plot(soln_1.x, soln_1.y, IEM_X, IEM_Y);
legend('ODE45', 'Improved Euler', 'Location', 'Best');
title("y' = 1/(y^2)");
ylabel('y');
xlabel('t');
% Major differences:
% At the beginning of the solutions, the IEM solver provides a smoother 
% solution, which is closer to the exact solution.
% However, the two solutions produced are relatively close.

% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
f_1 = @(t,y) 1-((t*y)/2);
t0 = 0;
tN = 10;
y0 = -1;
soln_1 = ode45(f_1,[t0, tN], y0);
[IEM_X, IEM_Y] = IEM_solver(f_1,t0,tN,y0,0.01);
subplot(2,2,3);
plot(soln_1.x, soln_1.y, IEM_X, IEM_Y);
legend('ODE45', 'Improved Euler', 'Location', 'Best');
title("y' = 1-((t*y)/2)");
ylabel('y');
xlabel('t');
% Major differences:
% The solution produced by the IEM solver is smoother than the one
% approximated by ODE45. This is appearant near t = [1,3]. This is likely
% due to the reason that the step size is smaller for the IEM solver.

% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
f_1 = @(t,y) (y^3)-(t^2);
t0 = 0;
tN = 10;
y0 = 1;
soln_1 = ode45(f_1,[t0, tN], y0);
[IEM_X, IEM_Y] = IEM_solver(f_1,t0,tN,y0,0.01);
subplot(2,2,4);
plot(soln_1.x, soln_1.y, IEM_X, IEM_Y);
legend('ODE45', 'Improved Euler', 'Location', 'Best');
title("y' = (y^3)-(t^2)");
ylabel('y');
xlabel('t');
% Major differences:
% Within the given time frame, it is hard to spot any difference. Both of
% the graphs seem to be appraching positive infinity after around t = 5.2.
% (However, ODE45 gives a warning stating that it cannot meet the
% integration tolerances without reducing the step size smaller than it is
% allowed around t = 5.07. This means that it is increasing too rapidly.)

% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
clear; close all; clc;
f = @(t,y) 2*t*sqrt(1-(y^2));
T = linspace(0,0.5);
y0 = 0;
EM_soln_y = euler(f,y0,T);
subplot(1,1,1);
plot(T,EM_soln_y);
title("Euler's Method: y' = 2*t*sqrt(1-(y^2))");
ylabel('y');
xlabel('t');
%
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
% The exact solution of the IVP:  y = sin(t^2)
fprintf('The analytical solution shows that y(0.5) = %f. \n', sin(0.5^2));
%
% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
% E_n = ((1+M)*dt/2)*(exp(M*dt)-1)*N
% M: Maximum of f, df/dt, df/dy
% The intervals for t and y are [0,0.5] and [0, 0.25] respectively
% M = 2 in this case, since the maximums of f, df/dt, df/dy all do not
% exceed 2
% dt: (time) step size, which is 0.005 in this case
% N: number of steps, which is 100 in this case
%
% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
N = 100;
M = 2;
dt = 0.005;
E_n = ((1+M)*dt/2)*(exp(M*dt)-1)*N;
fprintf('With a step size of 0.005, the error bound is %f.\n', E_n);
fprintf('The actual error is %f.\n', abs(sin(0.5^2) - EM_soln_y(N)));
%
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.
% Use 150 steps instead of 100, i.e. N = 150
% As a result, with t0 = 0 and tN = 0.5, dt = 0.0025
N = 200;
M = 2;
dt = 0.0025;
T = linspace(0,0.5,200);
EM_soln_y = euler(f, y0, T);
E_n = ((1+M)*dt/2)*(exp(M*dt)-1)*N;
fprintf('With a step size of 0.025, the error bound is %f.\n', E_n);
fprintf('The actual error is %f.\n', abs(sin(0.5^2) - EM_soln_y(N)));

% For large N values, the error is expected to be O(n)
% As we increase dt by 2, the error reduces by half the original, i.e. from 
% 0.007538 to 0.003759, which is what we expected.


%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

% The formula for updating the step size is attempting to make sure that
% with the step size it is using, the result will be within the tolerance 
% range, without having a step size too small. 0.9 is a reasonable ratio it
% reduces of the original step size. With the formula, if the error is much
% larger than the tolerance range, the step size will be greatly reduced. 
% On the other hand, if the error is close to the tolerance range, it will
% only change slightly.
%
%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
clear; close all; clc;
f = @(t,y) 2*t*sqrt(1-(y^2));
y0 = 0;
T = linspace(0,0.75,30);
h = 0.025;
EM_soln_y = euler(f,y0,T);
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
f = @(t,y) 2*t*sqrt(1-(y^2));
y0 = 0;
t0 = 0;
tN = 0.75;
h = 0.025;
[AEM_X, AEM_Y] = AEM_solver(f, t0, tN, y0, h);
%
% (c) Plot both approximations together with the exact solution.
Soln_exact_y = sin(T.^2);
Soln_exact_x = linspace(0,0.75,30);
subplot(1,1,1);
plot(T, EM_soln_y, AEM_X, AEM_Y, Soln_exact_x, Soln_exact_y,'x', 'MarkerSize',10, 'LineWidth', 2);
title("y' = 2*t*sqrt(1-(y^2))");
legend("Euler's Method", "Adaptive Euler's Method", "Exact Solution", 'Location', 'Best');
ylabel('y');
xlabel('t');

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.

% The Adpative Euler's Method solver (AEM) is closer to the exact solution. The reason
% being, the number of steps that it took was way larger than the Euler's
% Method (EM). While the AEM took 5600+ steps, the EM only took 30, which
% can be observed in the "Workspace".
 
% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
clear; close all; clc;
f = @(t,y) 2*t*sqrt(1-(y^2));
y0 = 0;
t0 = 0;
tN = 1.5;
h = 0.025;
T = linspace(0,1.5,30);
[AEM_X, AEM_Y] = AEM_solver(t0, tN, y0, h,f);
EM_soln_y = euler(f,y0,T);
[AEM_X, AEM_Y] = AEM_solver(f, t0, tN, y0, h);
Soln_exact_y = sin(T.^2);
Soln_exact_x = linspace(t0,tN,30);
plot(T, EM_soln_y, AEM_X, AEM_Y, Soln_exact_x, Soln_exact_y,'x', 'MarkerSize',10, 'LineWidth', 2);
title("y' = 2*t*sqrt(1-(y^2))");
legend("Euler's Method", "Adaptive Euler's Method", "Exact Solution", 'Location', 'Best');
ylabel('y');
xlabel('t');
%
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

% The exact solution and the Adaptive Euler's Method (AEM) and ODE45
% approximated solutions are very different starting from when t is around 
% pi/2. This is because of the term sqrt(1-y^2) in the equation. When the
% approximations provide a result where y is larger than 1, the square root
% term becomes imaginary. Since MATLAB only considers real numbers, it 
% contributes to the graphs looking odd and further away from the exact 
% solution.
