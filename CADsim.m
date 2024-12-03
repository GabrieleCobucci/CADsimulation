function v = CADsim(d,r,m,nUnit,F,rho)

%-------------------------------------------------------------------------%
%This function computes the r-dimensional classical simulation of a
%d-dimensional quantum ensemble of pure states subject to white noise,
%i.e. rho_x = v*psi_x + (1-v)/d*Id

%Inputs:
% - d: physical dimension of the states in the ensemble
% - r: dimensionality of the classical simulation
% - m: number of states in the ensemble
% - nUnit: number of unitaries to simulate the ensemble
% - F: unitaries to simulate the ensemble
% - rho: d-dimensional quantum ensemble {rho_x}

%Output:
% - v: maximum visibility for which an r-dimensional classical simulation
% is possible
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

%Number of simulation boxes
nSub = nchoosek(d,r); %Number of subspaces
subsets = nchoosek([1:d],r); %Subspaces

nLambda = nUnit*nSub; %Number of boxes

%Construction of projectors

%Subspace projectors
for y = 1 : nUnit
    for mu = 1 : nSub
        subarr = zeros(1,d);
        subarr(1,subsets(mu,:)) = 1; %Fixing to 1 the interested subset
        proj{y,mu} = F{y}*diag(subarr)*F{y}'; %Projectors
    end
end

%Vectors in the subspace bases
for y = 1 : nUnit
    for mu = 1 : nSub
        L = subsets(mu,:);
        for l = 1 : length(L)
            subvec{l,y,mu} = F{y}*id(:,L(l));
        end
    end
end

%% SDP

%Constraints
C = [];

%Variables
v = sdpvar(1); %Visibility

for x = 1 : m
    for y = 1 : nUnit
        for mu = 1 : nSub
            for k = 1 : r
                p{k,x,y,mu} = sdpvar(1); %Unnormalized probabilities in the diagonal decomposition of the states from the simulation boxes
                C = [C, p{k,x,y,mu} >= 0];
            end
        end
    end
end

%Simulation coefficients
sumq = 0;
for y = 1 : nUnit
    for mu = 1 : nSub
        q{y,mu} = sdpvar(1);
        C = [C, q{y,mu} >= 0];
        sumq = sumq + q{y,mu};
    end
end
C = [C, sumq == 1];

for x = 1 : m
    for y = 1 : nUnit
        for mu = 1 : nSub
            sump = 0;
            for k = 1 : r
                sump = sump + p{k,x,y,mu};
            end
            C = [C, sump == q{y,mu}]; %Normalization constraint
        end
    end
end

%Explicit definition of the states from the simulation boxes
for x = 1 : m
    tau{x} = 0;
    for y = 1 : nUnit
        for mu = 1 : nSub
            for k = 1 : r
                tau{x} = tau{x} + p{k,x,y,mu}*subvec{k,y,mu}*subvec{k,y,mu}';
            end
        end
    end
end

%White noise
for x = 1 : m
    C = [C, v*rho{x} + (1-v)/d*id == tau{x}];
end


%SolveSDP
disp('Options')
ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
diagnostic=solvesdp(C,-v,ops)

v=double(v)

end