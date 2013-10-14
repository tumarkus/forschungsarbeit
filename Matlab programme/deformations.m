





function [ deform ] = deformations(accj, forcel, forcer, PROP1, PROP2, PROP3, PROP4, PAR, teta, taul, taur)
%deformations Calculates the link deformations from a given angular acceler
%   Detailed explanation goes here
forcejlx=forcel(1);
forcejly=forcel(2);
for i=1:size(PROP2,1)
    incx=PROP2(i,2)*PAR(2,2)*(accj(1,1)+(accj(2,1)-accj(1,1))*PROP2(i,1));
    incy=PROP2(i,2)*PAR(2,2)*(accj(1,2)+(accj(2,2)-accj(1,2))*PROP2(i,1));
    forcejlx=forcejlx+incx;
    forcejly=forcejly+incy;
end
forcejrx=forcer(1);
forcejry=forcer(2);
for i=1:size(PROP4,1)
    incx=PROP4(i,2)*PAR(4,2)*(accj(3,1)+(accj(4,1)-accj(3,1))*PROP4(i,1));
    incy=PROP4(i,2)*PAR(4,2)*(accj(3,2)+(accj(4,2)-accj(3,2))*PROP4(i,1));
    forcejrx=forcejrx+incx;
    forcejry=forcejry+incy;
end
forcejlx
forcejly
%% calculating K1-K4
K1=zeros(2*size(PROP1,1));
K2=zeros(2*size(PROP2,1));
K3=zeros(2*size(PROP3,1));
K4=zeros(2*size(PROP4,1));

for i=1:size(PROP1,1)-1
    l=PAR(1,1)*(PROP1(i+1,1)-PROP1(i,1));
    k=PROP1(i,3)*PROP1(i,4)/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
    k_i=zeros(2*size(PROP1,1));
    k_i(2*i-1 : 2*i+2, 2*i-1 : 2*i+2)=k;
    K1=K1+k_i;
end

for i=1:size(PROP2,1)-1
    l=PAR(2,1)*(PROP2(i+1,1)-PROP2(i,1));
    k=PROP2(i,3)*PROP2(i,4)/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
    k_i=zeros(2*size(PROP2,1));
    k_i(2*i-1 : 2*i+2, 2*i-1 : 2*i+2)=k;
    K2=K2+k_i;
end

for i=1:size(PROP3,1)-1
    l=PAR(3,1)*(PROP3(i+1,1)-PROP3(i,1));
    k=PROP3(i,3)*PROP3(i,4)/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
    k_i=zeros(2*size(PROP3,1));
    k_i(2*i-1 : 2*i+2, 2*i-1 : 2*i+2)=k;
    K3=K3+k_i;
end

for i=1:size(PROP4,1)-1
    l=PAR(4,1)*(PROP4(i+1,1)-PROP4(i,1));
    k=PROP4(i,3)*PROP4(i,4)/l^3*[12, 6*l, -12, 6*l; 6*l, 4*l^2, -6*l, 2*l^2; -12, -6*l, 12, -6*l; 6*l, 2*l^2, -6*l, 4*l^2];
    k_i=zeros(2*size(PROP4,1));
    k_i(2*i-1 : 2*i+2, 2*i-1 : 2*i+2)=k;
    K4=K4+k_i;
end

%% Calculating F1vek to F4vek
F1vek=zeros(2*size(PROP1,1),1);
F2vek=zeros(2*size(PROP2,1),1);
F3vek=zeros(2*size(PROP3,1),1);
F4vek=zeros(2*size(PROP4,1),1);
for i=1:(size(PROP1,1))
   Finc=zeros(2*size(PROP1,1),1);
   f=[-PAR(1,2)*PROP1(i,2)*norm([sin(teta(1)); -cos(teta(1))]*(accj(1,:)*PROP1(i,1))); 0];
   Finc(2*i-1:2*i)=f;
   F1vek=F1vek+Finc;
end
F1vek(2*size(PROP1,1)-1: 2*size(PROP1,1))=F1vek(2*size(PROP1,1)-1: 2*size(PROP1,1))+[-[sin(teta(1)), -cos(teta(1))]*[forcejlx; forcejly]; taul(2)];
F1vek
for i=1:(size(PROP2,1))
   Finc=zeros(2*size(PROP2,1),1);
   f=[-PAR(2,2)*PROP2(i,2)*norm([sin(teta(1)+teta(2)); -cos(teta(1)+teta(2))]*((accj(1,:)+PROP2(i,1)*((accj(2,:)-accj(1,:)))*PROP2(i,1)))); 0];
   Finc(2*i-1:2*i)=f;
   F2vek=F2vek+Finc;
end

F2vek(2*size(PROP2,1)-1: 2*size(PROP2,1))=F2vek(2*size(PROP2,1)-1: 2*size(PROP2,1))+[-[sin(teta(1)+teta(2)), -cos(teta(1)+teta(2))]*forcel'; 0];
F2vek
for i=1:(size(PROP3,1))
   Finc=zeros(2*size(PROP3,1),1);
   f=[-PAR(3,2)*PROP3(i,2)*norm([sin(teta(3)); -cos(teta(3))]*(accj(3,:)*PROP3(i,1))); 0];
   Finc(2*i-1:2*i)=f;
   F3vek=F3vek+Finc;
end
F3vek(2*size(PROP3,1)-1: 2*size(PROP3,1))=F3vek(2*size(PROP3,1)-1: 2*size(PROP3,1))+[-[sin(teta(3)), -cos(teta(3))]*[forcejrx; forcejry]; taur(2)];
F3vek
for i=1:(size(PROP4,1))
   Finc=zeros(2*size(PROP4,1),1);
   f=[-PAR(4,2)*PROP4(i,2)*norm([sin(teta(3)+teta(4)); -cos(teta(3)+teta(4))]*((accj(4,:)+PROP4(i,1)*((accj(4,:)-accj(3,:)))*PROP4(i,1)))); 0];
   Finc(2*i-1:2*i)=f;
   F4vek=F4vek+Finc;
end
F4vek(2*size(PROP4,1)-1: 2*size(PROP4,1))=F4vek(2*size(PROP4,1)-1: 2*size(PROP4,1))+[-[sin(teta(3)+teta(4)), -cos(teta(3)+teta(4))]*[forcer']; 0];
F4vek
%% Calculation of the deformations
%ATTENTION: Moments are in wrong x-y direction!! ToDo: FIX
%ToDo:Check if Moments at encastre are same as calculated

V1vek=zeros(2*size(PROP1,1),1);
V1vek(3:2*size(PROP1,1))=K1(3:2*size(PROP1,1),3:2*size(PROP1,1))^-1*F1vek(3:2*size(PROP1,1))
V2vek=zeros(2*size(PROP2,1),1);
V2vek(3:2*size(PROP2,1))=K2(3:2*size(PROP2,1),3:2*size(PROP2,1))^-1*F2vek(3:2*size(PROP2,1))
V3vek=zeros(2*size(PROP3,1),1);
V3vek(3:2*size(PROP3,1))=K3(3:2*size(PROP3,1),3:2*size(PROP3,1))^-1*F3vek(3:2*size(PROP3,1))
V4vek=zeros(2*size(PROP4,1),1);
V4vek(3:2*size(PROP4,1))=K4(3:2*size(PROP4,1),3:2*size(PROP4,1))^-1*F4vek(3:2*size(PROP4,1))
size(V1vek)
deform=[V1vek(size(V1vek,1)-1), V1vek(size(V1vek,1));-V2vek(size(V2vek,1)-1), -V2vek(size(V2vek,1));V3vek(size(V3vek,1)-1), V3vek(size(V3vek,1));-V4vek(size(V4vek,1)-1), -V4vek(size(V4vek,1))]

end
