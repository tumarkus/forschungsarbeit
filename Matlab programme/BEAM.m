function [ erg] = BEAM( number, E, I, L, q0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
l=L/number;
k=E*I/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
K=zeros (2*(number+1));

for i=0:number-1
    
    k_i=zeros(2*(number+1));
    k_i(2*i+1 : 2*i+4, 2*i+1 : 2*i+4)=k;
    K=K+k_i;

end
F_vek=zeros(2*(number+1), 1);
f=[q0*l/2; 0; q0*l/2;0];
for i=0:number-1
    f_i=zeros(2*(number+1), 1);
    f_i(2*i+1:2*i+4, 1)=f;
    F_vek=F_vek+f_i;
end
v_vek=zeros(2*(number+1), 1);
v_vek(3:(2*(number+1)))=K(3:(2*(number+1)), 3: (2*(number+1)))^-1*F_vek(3:2*(number+1));
x=0:l:L;
x=x';
y=zeros(number+1, 1);
for i=1:number+1
    y(i)=v_vek(2*i-1);
end


erg=[x, y];
