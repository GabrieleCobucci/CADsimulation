function vis = unitarymax(x)

%Ensemble parameters
m = 7; %Number of states
d = 7; %Physical dimension
r = 7; %CAD

%Number of unitaries for simulation
nUnit = 20;

%Physical identity
id = eye(d);

%Unitaries
for y = 1 : nUnit
    lambda = zeros(d);
    for i = 1 : d
        for j = 1 : d
            lambda(i,j) = x(d*(y-1)+i*d-d+j);
        end
    end
    F{y} = UC(lambda);
end

%% Ensemble construction
for l = 1 : d-1
    rho{l} = id(:,l)*id(:,l)';
end

four = 0;
for l = 0 : d-1
    four = four + 1/sqrt(d)*id(:,l+1);
end

rho{d} = four*four';

%% SDP
v = CADsim(d,r,m,nUnit,F,rho);

%Maximizing with fmincon
vis = -v

end