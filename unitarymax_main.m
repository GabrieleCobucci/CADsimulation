clear all
clc

%Dimension
d = 5;

%Constraints
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
lb = [];
ub = [];

dat = 0; %Results for each iteration
itnum = 50; %Number of iterations

for num = 1 : itnum
    %Random unitary to start the optimization
    x0 = rand(d^2,1);
    
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp' );
    [x,unitmax] = fmincon(@unitarymax,x0,A,b,Aeq,beq,lb,ub);

    dat = [dat, unitmax];
end

%Minimum over iterations
min(dat)