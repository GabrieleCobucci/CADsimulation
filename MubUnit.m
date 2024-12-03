function U = MubUnit(d)

%-------------------------------------------------------------------------%
%This function generates the set of unitaries from the computational basis
%to the d+1 MUBs in a d-dimensional Hilbert space, when d is odd.

%Input:
% - d: dimension of the Hilbert space

%Output:
% - U: set of unitaries, U{1},..., U{d+1};
%-------------------------------------------------------------------------%

%Physical identity
id = eye(d);
omega = exp(2*pi*i/d);

%MUBs
for l = 0 : d-1
    mub{1,l+1} = id(:,l+1);
end

for j = 2 : d+1
    for l = 0 : d-1
        mub{j,l+1} = 0;
        for k = 0 : d-1
            mub{j,l+1} = mub{j,l+1} + 1/sqrt(d)*omega^(k*l + (j-2)*k^2)*id(:,k+1);
        end
    end
end

for y = 1 : d+1
    U{y} = 0;
    for k = 0 : d-1
        U{y} = U{y} + mub{y,k+1}*id(:,k+1)';
    end
end

end