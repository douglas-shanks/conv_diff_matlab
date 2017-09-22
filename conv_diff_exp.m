function [A]=conv_diff_exp(N,R,P)
%% conv_diff_exp.m
%   u_{t} - \alpha div u  - \beta \Delta u = 0, in Omega = [0,1]^2
%                                   u(x,0) = u_{0}
%                                        u = 0 on \partial\Omega
% 
dx=(1)./(N-1);
dx2=dx^2;
a=zeros(N*7);   % x's
b=zeros(N*7);   % y's
c=zeros(N*7);   % nz's
pos=1;

%% form A
for m = 1: N
    for j = 1: N
        
        %% Interior
            a(pos) = (m-1)*(N)+j;
            b(pos) = (m-1)*(N)+j;
            c(pos)=1.0 - 4.0*R;
            pos = pos+1;
        
        %% NSEW in stencil
        if j<N
            a(pos) = (m-1)*(N)+j;
            b(pos) = (m-1)*(N)+j+1;
            c(pos)=(R + R*P); %% N
            pos = pos+1;
        end
        if j>1
            a(pos) = (m-1)*(N)+j;
            b(pos) = (m-1)*(N)+j-1;
            c(pos)=(R+ R*P); %% S
            pos = pos+1;
        end
        if m<N
            a(pos) = (m-1)*(N)+j;
            b(pos) = m*(N)+j;
            c(pos)=(R - R*P); %% E
            pos = pos+1;
        end
        if m>1
            a(pos) = (m-1)*(N)+j;
            b(pos) = (m-2)*(N)+j;
            c(pos)=(R - R*P); %% W
            pos = pos+1;
        end

    end
a=a(1:pos-1);
b=b(1:pos-1);
c=c(1:pos-1);
A=sparse(a,b,c);

end
