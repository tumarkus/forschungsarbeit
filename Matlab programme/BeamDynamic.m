function [ output_args ] = BeamDynamic( E,I,L,number,Mass,Fm,freq,timeinc,endtime,sloMo)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%Constants:
g=9.81;
l=L/number;
%Variables:
time=0;
Q=zeros(2*(number+1),1);
Qdot=zeros(2*(number+1),1);
Qddot=zeros(2*(number+1),1);
V=zeros(number+1,1);

%K-Matrix as always
k=E*I/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
K=zeros (2*(number+1));
for i=0:number-1
    
    k_i=zeros(2*(number+1));
    k_i(2*i+1 : 2*i+4, 2*i+1 : 2*i+4)=k;
    K=K+k_i;

end
%M is the plumbed Mass-Matrix describing how much mass is attached to every
%node
m=Mass/2/number*[1,0,0,0;0,0,0,0;0,0,1,0;0,0,0,0];
M=zeros(2*(number+1));
for i=0:number-1
    
    m_i=zeros(2*(number+1));
    m_i(2*i+1 : 2*i+4, 2*i+1 : 2*i+4)=m;
    M=M+m_i;

end
%the C-Matrix describes the damping. 0.1 is just a guess.
C=zeros(2*(number+1));
for i=1:((number+1))
    C(2*i-1, 2*i-1)=20;
end
%F_vek is the force-vector caused by gravitation
F_vek=zeros(2*(number+1), 1);
for i=1:(2*(number+1))
    F_vek(i)=M(i,i)*g;
end
%Kmod is a modified K-Matrix used for calculating Q-vektor
Kmod=(4/timeinc^2*M+2/timeinc*C+K);
%Here starts Neumanns-Algorithm
while (time<endtime)
    %Here I modify the f-vektor with a sinus force on the right boundary
    F_vek2=F_vek;
    F_vek2(end-1)=F_vek2(end-1)+sin(2*3.1415*freq*time)*Fm;
    Fmod=F_vek2+M*(4/timeinc^2*Q+4/timeinc*Qdot+Qddot)+C*(2/timeinc*Q+Qdot);
    Qold=Q;
    Qdotold=Qdot;
    Qddotold=Qddot;
    Q(3:end)=Kmod(3:end,3:end)\Fmod(3:end);
    Q(1:2)=[0;0];
    Qddot=4/timeinc^2*Q-4/timeinc^2*Qold-4/timeinc*Qdotold-Qddotold;
    Qdot=Qdotold+timeinc/2*Qddotold+timeinc/2*Qddot;
    
    %calculating the displacements
    for i=1:(number+1)
        V(i)=Q(2*i-1);
    end
    
    plot((1:(number+1)),V)
    axis ([0 number+1 -0.005 0.01])
    pause(timeinc*sloMo);
    clc;
    time=time+timeinc
end


end

