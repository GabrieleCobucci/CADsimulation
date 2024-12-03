%-------------------------------------------------------------------------%
%Example of classical simulation search by using random unitaries.
%-------------------------------------------------------------------------%

clear all
clc

res = 0; %Results from each iterations
itnum = 50; %Number of iterations

for num = 1 : itnum

%Ensemble parameters
m = 5; %Number of states
d = 5; %Physical dimension
r = 5; %QAD

%Number of unitaries for simulation
nUnit = 20;

%Physical identity
id = eye(d);

%Random unitaries
for y = 1 : nUnit
    U{y} = RandomUnitary(d);
end

%Ensemble
for x = 1 : d-1
    rho{x} = id(:,x)*id(:,x)';
end

s = 0;
for k = 1 : d
    s = s + 1/sqrt(d)*id(:,k);
end

rho{d} = s*s';

%Max visibility
v = CADsim(d,r,m,nUnit,U,rho);

res = [res, v];

end