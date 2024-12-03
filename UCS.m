%
%   COMPOSITE PARAMETERIZATION of SU(d) /// version 1.5 /// 12.04.2011
%
%   (c) Christoph Spengler 2011, Faculty of Physics, University of Vienna
%       Contact: christoph.spengler@univie.ac.at
%
%   Usage : UCS(lambda)
%
%   lambda - dxd real matrix (lambda(d,d) is ignored)
%   lambda(a,b) diagonal components a=b with a,b in {1,..,d-1} - absolute phases for a in [0,2*pi]
%   lambda(a,b) upper right components a<b - rotations in a-b plane in [0,pi/2]
%   lambda(a,b) lower left components a>b - relative phases between a-b in [0,pi]
%
%   References: --- PLEASE CITE THESE PAPERS WHEN USING THIS FILE ---
%
%   Ch.Spengler, M.Huber, B.C.Hiesmayr
%   'A composite parameterization of unitary groups, density matrices and subspaces'
%   arXiv:1004.5252 // J. Phys. A: Math. Theor. 43, 385306 (2010) 
%
%   Ch.Spengler, M.Huber, B.C.Hiesmayr
%   'Composite parameterization and Haar measure for all unitary and special unitary groups'
%   arXiv:1103.3408
%

function unitary=UCS(lambda)

d=length(lambda);

unitary=1;
ex2=1;

for m=(d-1):-1:1
    ex1(1,d-m)=0;
    ex2(d-m+1,1)=0;
    unitary=[ex1;unitary];
    unitary=[ex2 unitary];    
    for n=(d-m+1):-1:2 
        A=eye(d-m+1);
        A(1,1)=exp(1i*lambda(n+m-1,m))*cos(lambda(m,n+m-1));
        A(n,n)=exp(-1i*lambda(n+m-1,m))*cos(lambda(m,n+m-1));
        A(n,1)=-exp(-1i*lambda(n+m-1,m))*sin(lambda(m,n+m-1));
        A(1,n)=exp(1i*lambda(n+m-1,m))*sin(lambda(m,n+m-1));
        unitary=A*unitary;
    end
end

for k=1:d-1
    unitary(:,k)=unitary(:,k)*exp(1i*lambda(k,k));
    unitary(:,d)=unitary(:,d)*exp(-1i*lambda(k,k));
end


end
