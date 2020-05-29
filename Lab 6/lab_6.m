%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|. Also in this lab, you will write your own
% ODE solver using Laplace transforms and check whether the result yields
% the correct answer.
%
% You will learn how to use the |laplace| routine. 
% 
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in the template, including appropriate descriptions 
% in each step. Save the m-file and submit it on Quercus.
%
% Include your name and student number in the submitted file.
%
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, based on 
% MAT292, Fall 2013, Sinnamon & Sousa

%% Using symbolic variables to define functions
% 
% Recall the use of symbolic variables and function explained in the MATLAB
% assignment #2.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.

% Clearing all previous work
clc; close all; clear;

% Definition of symbolic variables
syms t s a;

% Definition of the function
f(t) = exp(2*t)*t^3;

% Computation of the Laplace Transform of the function
F = laplace(f);

% Displaying the Laplace Transform of the function
disp(F);

% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|

% Definition of the function
g = ((s - 1)*(s - 2))/(s*(s + 2)*(s - 3));

% Computation of the Laplace Transform of the function
G = ilaplace(g);

% Displaying the Laplace Transform of the function
disp(G);

% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 

% Definition of symbolic variables
syms a f(t) g(t) t

% Definition of the new function
g(t) = exp(a*t)*f(t);

% Computation of the Laplace Transform of the functions
F = laplace(f(t));
G = laplace(g(t));

% Displaying the Laplace Transform of the functions
disp(F);
disp(G);

% Explaination: 
% f(t) is an arbitrarily function of t which has the Laplace Transform of 
% F(s) = laplace(f(t), t, s). g(t) is exp(at)*f(t), and its Laplace 
% Transform is G = laplace(f(t), t, s - a). By comparison, we can
% see that G = F(s-a) since in the last input of laplace(f(t), t, s), s is 
% replaced by (s-a) in G = laplace(f(t), t, s - a). 

% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.


%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

% Clearing all previous work
clc; close all; clear;

% Definition of symbolic variables
syms t s a;

% Assigning a value to a
a = 3;

% Definition of Heaviside function |u_a(t)| at |a|
u_a(t) = heaviside(t-a);

% Definition of the functions
f(t) = exp(2*t)*t^3;
g(t) = u_a(t)*f(t-a);

% Computation of the Laplace Transform of the functions
F = laplace(f(t));
G = laplace(g(t));

% Displaying the Laplace Transform of the functions
disp(F);
disp(G);

% The Laplace Transform of f(t) is F(s) = 6/(s - 2)^4, and the Laplace
% Transform of g(t) = u_a(t)*f(t-a) is G(s) = (6*exp(-3*s))/(s - 2)^4, 
% which is exp(-3*s)*F(s). In the textbook, Theorem 5.5.1 states that the
% Laplace Transform of u_c(t)*f(t-c) is equal to exp(-c*t) multiplied by
% the Laplace Transform of f(t), which is F(s). In this case, c = a = 3.
% Therefore, the result is confirmed by the theorem.

%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

% Clearing all previous work
clc; close all; clear;

% Definition of the unknown function, its variable and the Laplace 
% tranform of the unknown
syms y(t) t Y s;

% Definition of the ODE
ODE = diff(y(t),t,3)+2*diff(y(t),t,2)+diff(y(t),t,1)+2*y(t)+cos(t) == 0;

% Computation of the Laplace transform of the ODE
L_ODE = laplace(ODE);

% Use/Sub in the initial conditions
L_ODE = subs(L_ODE,y(0),0);
L_ODE = subs(L_ODE,subs(diff(y(t), t), t, 0),0);
L_ODE = subs(L_ODE,subs(diff(y(t), t, 2), t, 0),0);

% Factor out the Laplace transform of |y(t)|
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y = solve(L_ODE,Y);

% Get the inverse Laplace transform (which is the solution to the original
% IVP)
y = ilaplace(Y);

% Plot the solution
ezplot(y,[0,10*pi]);

% Label the graph
title("Solution of y'''+2y''+y'+2*y=-cos(t)");
ylabel('y');
xlabel('t');

% There is not an initial condition for which |y| remains bounded as |t| 
% goes to infinity

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |14|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

% Clearing all previous work
clc; close all; clear;

% Definition of the unknown function, its variable and the Laplace 
% tranform of the unknown
syms y(t) t Y s;

% Definition of Heaviside function |u_0(t)| at |0|/|u_2(t)| at |2|
% and |u_5(t)| at |5|
u_0(t) = heaviside(t);
u_2(t) = heaviside(t-2);
u_5(t) = heaviside(t-5);

% Definition of the piecewise function g(t)
g(t) = 3*u_0(t)+(t-2)*u_2(t)+(-t+4)*u_5(t);

% Definition of the ODE
ODE = diff(y(t),t,2)+2*diff(y(t),t)+5*y(t)-g(t) == 0;

% Computation of the Laplace transform of the ODE
L_ODE = laplace(ODE);

% Use/Sub in the initial conditions
L_ODE=subs(L_ODE,y(0),2);
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),1);

% Factor out the Laplace transform of |y(t)|
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y=solve(L_ODE,Y);

% Get the inverse Laplace transform (which is the solution to the original
% IVP)
y = ilaplace(Y);

% Plot the solution
ezplot(y,[0,12,0,2.25]);

% Label the graph
title("Solution of y''+2y'+5y = g(t)");
ylabel('y');
xlabel('t');

%% Exercise 5
%
% Objective: Solve an IVP with a periodic function using the Laplace
% transform.
%
% Details:
%
% * Follow the instructions on the website
%  http://instruct.math.lsa.umich.edu/lecturedemos/ma216/docs/7_5/
%
% * Do the part labelled |Outside of Lecture|
%
% * Note that |u(t-a)| is the Heaviside function |u_a(t)| defined 
% in our textbook
%
% * Check Theorem 5.5.3 (page 349) to know how to define the Laplace 
% transform of a periodic function like the one in this exercise (and check
% the function |int| on MATLAB for symbolic integration).
%
% * Hint: Use only the second-order DE given on the linked website
% * Hint: You should use the numbers provided (the less precise ones)
% * Hint: You should obtain an exact solution. In the process, you should need symbolic integration(s).

% Clearing all previous work
clc; close all; clear;

% Definition of the unknown function, its variable and the Laplace 
% tranform of the unknown
syms y(t) Y f(t) t x s

% Definition of all given parameters (from the website)
k1 = 0.12;
k2 = 0.38;
k3 = 0.04;
I0 = 24;        % magnitude of dosing
t0 = 6;         % length of dosing
t1 = 24;        % time between doses
J = I0/t0;

% Definition of Heaviside function |u_0(t)| at |0|/|u_6(t)| at |6|
u_0(t) = heaviside(t);
u_6(t) = heaviside(t-6);

% Computation of I
f(t) = J*(u_0(t)-u_6(t));
INT = int(f(t)*exp(-s*t),t,0,t1);

% Getting the Laplace Transform of I
I_Laplace = INT/(1-exp(-s*t1));

% Definition of the ODE of the left hand side of the equation
ODE_L = diff(y(t),t,2)+(k1+k2+k3)*diff(y(t),t)+(k3*(k1+k2))*y(t);

% Computation of the Laplace Transform of the ODE of the left hand side 
% of the equation
L_ODE_Laplace = laplace(ODE_L);

% Equating the left hand side and the right hand side of the equation
L_ODE = L_ODE_Laplace == k1*I_Laplace;

% Use/Sub in the initial conditions
L_ODE=subs(L_ODE,y(0),0);
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),0);

% Factor out the Laplace transform of |y(t)|
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y=solve(L_ODE,Y);

% Get the inverse Laplace transform (which is the solution to the original
% IVP)
y = ilaplace(Y);

% Plot the solution
ezplot(y,[0,250]);

% Label the graph
title("Solution");
ylabel('y');
xlabel('t');
