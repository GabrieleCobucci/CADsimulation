clear all
clc

%Dimension
d = 7;

%Constraints
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
lb = [];
ub = [];

res = 0; %Results from each iterations
itnum = 50; %Number of iterations

%Number of unitaries
nUnit = 20;

for num = 1 : itnum
    %Random unitary to start the optimization
    x0 = rand(nUnit*d^2,1);

    options = optimoptions('fmincon','Display','iter','Algorithm','sqp' );
    [x,unitmax] = fmincon(@unitarymax,x0,A,b,Aeq,beq,lb,ub);

    res = [res, unitmax];
end

%Minimum over iterations
min(res)

