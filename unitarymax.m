function vis = unitarymax(x)

%Ensemble parameters
m = 3; %Number of states
d = 5; %Physical dimension
r = 5; %CAD

%Number of unitaries for simulation
nUnit = d+1;

%Number of MUBs
nMub = d+1;

%Physical identity
id = eye(d);

%MUBs unitaries
U = MubUnit(d);

%Rotation construction
lambda = zeros(d);
for i = 1 : d
    for j = 1 : d
        lambda(i,j) = x(i*d-d+j);
    end
end
A = UC(lambda);

for y = 1 : nMub
    F{y} = A*U{y}*A';
end


%% Ensemble construction

for l = 0 : d-1
    rho{l+1} = id(:,l+1)*id(:,l+1)';
end

four = 0;
for l = 0 : d-1
    four = four + 1/sqrt(d)*id(:,l+1);
end

rho{d} = four*four';

%% SDP

v = CADsim(d,r,m,nUnit,F,rho);

%FMINCON objective
vis = -v

end