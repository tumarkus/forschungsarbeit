function [taul, taur] = twoTimesTwoLink(l0, x, y, alpha, xend, yend, alphaend, time, PROP1,PROP2,PROP3,PROP4,PROPPL)
%twoTimesTwoLink.m
%calculates the trajectories of theta (tetatraj) of two given 2R
%manipulators with PROP1(end,3) and PROP2(end,3) for the arm length of the left arm and PROP3(end,3) and
%PROP4(end,3) as the arm length of the right arm for a given movement of the payload
%with lenght a from pose (x,y,alpha) to (xend, yend, alphaend).
%
%It also plots the trajectories of the two endeffectors.
PAR=[parameters(PROP1);parameters(PROP2);parameters(PROP3);parameters(PROP4);parameters(PROPPL)]
g=9.81;
steps=200;
timeinc=time/steps;
deform=[0,0;0,0;0,0;0,0];
%f1 and f2 are set to 1 if the attachment point of the payload is out of
%reach
f1=0; f2=0;
%I devide the differences in steps^2/2 increments, so that the position can
%be a quadratic function, the velocity linear and the accelerations
%constant
forcel=zeros(steps+1,2);
forcer=zeros(steps+1,2);

taul=zeros(steps+1, 2);
taur=zeros(steps+1, 2);


xinc=(xend-x)/(steps^2/2);
yinc=(yend-y)/(steps^2/2);
alphainc=(alphaend-alpha)/(steps^2/2);

v1=0; v2=0; v3=0; v4=0;

%xtraj is a vektor with the trajectory of the center of the payload


xtraj=zeros(steps+1, 1);
for i=1:steps/2+1
   xtraj(i)=x+xinc*(i-1)^2; 
end
for i=1:steps/2
    xtraj(steps+2-i)=xend-xinc*(i-1)^2;
end
%same with y
ytraj=zeros(steps+1, 1);
for i=1:steps/2+1
   ytraj(i)=y+yinc*(i-1)^2; 
end
for i=1:steps/2
    ytraj(steps+2-i)=yend-yinc*(i-1)^2;
end
%same with alpha
alphatraj=zeros(steps+1,1);
for i=1:steps/2+1
   alphatraj(i)=alpha+alphainc*(i-1)^2; 
end
for i=1:steps/2
    alphatraj(steps+2-i)=alphaend-alphainc*(i-1)^2;
end
%calculation of xdot, ydot and alphadot on every point of the trajectory
xdottraj=zeros(steps+1, 1);
for i=1:steps
    xdottraj(i)=(xtraj(i+1)-xtraj(i))/timeinc;
end
xdottraj(steps+1)=0;

ydottraj=zeros(steps+1, 1);
for i=1:steps
    ydottraj(i)=(ytraj(i+1)-ytraj(i))/timeinc;
end
ydottraj(steps+1)=0;

alphadottraj=zeros(steps+1, 1);
for i=1:steps
    alphadottraj(i)=(alphatraj(i+1)-alphatraj(i))/timeinc;
end
alphadottraj(steps+1)=0;
%same for accelerations
xddottraj=zeros(steps+1, 1);
for i=1:steps
    xddottraj(i)=(xdottraj(i+1)-xdottraj(i))/timeinc;
end
xddottraj(steps+1)=0;

yddottraj=zeros(steps+1,1);
for i=1:steps
    yddottraj(i)=(ydottraj(i+1)-ydottraj(i))/timeinc;
end
yddottraj(steps+1)=0;

alphaddottraj=zeros(steps+1, 1);
for i=1:steps
    alphaddottraj(i)=(alphadottraj(i+1)-alphadottraj(i))/timeinc;
end
alphaddottraj(steps+1)=0;

%putting position, velocity and acceleration in one matrix each
postraj=[xtraj, ytraj, alphatraj];
veltraj=[xdottraj, ydottraj, alphadottraj];
acctraj=[xddottraj, yddottraj, alphaddottraj];
%x1traj and y1traj are vektors for the trajectory of the left attachment
%point
x1traj=xtraj-cos(alphatraj)*PROPPL(end,3)*(PAR(5,1))
y1traj=ytraj-sin(alphatraj)*PROPPL(end,3)*(PAR(5,1))
%this for-loop checks if the left attachment point is out of reach
for n=1:1:steps+1
    if (x1traj(n)^2+y1traj(n)^2>=(PROP1(end,3)+PROP2(end,3))^2) f1=1;
    end;
end
%same with x2traj and y2traj
x2traj=xtraj+cos(alphatraj)*PROPPL(end,3)*(1-PAR(5,1))
y2traj=ytraj+sin(alphatraj)*PROPPL(end,3)*(1-PAR(5,1))
for n=1:1:steps+1
    if (((l0-x2traj(n))^2+y2traj(n)^2)>=(PROP1(end,3)+PROP2(end,3))^2) f2=1;
    end;
end
f1
f2
%here starts the calculation of the initial thetas
if ((f1==0) && (f2==0))
x1=x-cos(alpha)*PROPPL(end,3)*(PAR(5,1))
y1=y-sin(alpha)*PROPPL(end,3)*(PAR(5,1))
x2=x+cos(alpha)*PROPPL(end,3)*(1-(PAR(5,1)))
y2=y+sin(alpha)*PROPPL(end,3)*(1-(PAR(5,1)))
teta2=acos(-(PROP1(end,3)^2+PROP2(end,3)^2-x1^2-y1^2)/(2*PROP1(end,3)*PROP2(end,3)));
gamma1=atan2(y1,x1);
gamma2=acos((x1^2+y1^2+PROP1(end,3)^2-PROP2(end,3)^2)/(2*sqrt(x1^2+y1^2)*PROP1(end,3)));
teta1=gamma1-gamma2;
%
teta4=acos(-(PROP3(end,3)^2+PROP4(end,3)^2-(l0-x2)^2-y2^2)/(2*PROP3(end,3)*PROP4(end,3)));
gamma3=atan2((y2),(l0-x2));
gamma4=acos(((l0-x2)^2+y2^2+PROP3(end,3)^2-PROP4(end,3)^2)/(2*sqrt((l0-x2)^2+y2^2)*PROP3(end,3)));
teta3=3.1415926535898-gamma3-gamma4;
teta=[teta1; teta2; teta3; teta4];

%calculation of initial deformation
Jl=[-sin(teta(1))*PROP1(end,3)-sin(teta(1)+teta(2))*PROP2(end,3), -sin(teta(1)+teta(2))*PROP2(end,3); cos(teta(1))*PROP1(end,3)+cos(teta(1)+teta(2))*PROP2(end,3), cos(teta(1)+teta(2))*PROP2(end,3)];
Jr=[-sin(teta(3))*PROP3(end,3)-sin(teta(3)+teta(4))*PROP4(end,3), -sin(teta(3)+teta(4))*PROP4(end,3); cos(teta(3))*PROP3(end,3)+cos(teta(3)+teta(4))*PROP4(end,3), cos(teta(3)+teta(4))*PROP4(end,3)];
forcel0=[0,-PROPPL(end,4)*g/2];
forcer0=[0,-PROPPL(end,4)*g/2];
taul0=(Jl'*-(forcel(n,:))');
taur0=(Jr'*-(forcer(n,:))');
deform=deformations(zeros(4,2), forcel0, forcer0, PROP1, PROP2, PROP3, PROP4, PAR, teta, taul(n,:), taur(n,:));



%pos2 and pos4 are trajectories of the left and right endeffector, tetatraj
%is the trajectory of the theta vektor
pos2=zeros(steps+1, 2);
pos4=zeros(steps+1,2);
tetatraj=zeros(steps+1, 4);
linkveltraj=zeros(steps+1, 8);
%
n=0;
tetadot=[0;0;0;0];
while n<steps
    %here we calculate the position of the left and right endeffektor and save
    %it in pos2 and pos4
    n=n+1;
    pos2(n,1)=cos(teta(1))*PROP1(end,3)-sin(teta(1))*deform(1,1)+cos(teta(1)+teta(2)+deform(1,2))*PROP2(end,3)-sin(teta(1)+teta(2))*deform(2,1);
    pos2(n,2)=sin(teta(1))*PROP1(end,3)+cos(teta(1))*deform(1,1)+sin(teta(1)+teta(2)+deform(1,2))*PROP2(end,3)+cos(teta(1)+teta(2))*deform(2,1);
    pos4(n,1)=l0+cos(teta(3))*PROP3(end,3)-sin(teta(3))*deform(3,1)+cos(teta(3)+teta(4)+deform(3,2))*PROP4(end,3)-sin(teta(3)+teta(4))*deform(4,1);
    pos4(n,2)=sin(teta(3))*PROP3(end,3)+cos(teta(3))*deform(3,1)+sin(teta(3)+teta(4)+deform(3,2))*PROP4(end,3)+cos(teta(3)+teta(4))*deform(4,1);
    tetatraj(n,:)=teta';
    %x1dot and y1dot are the differces between the estimated position of the
    %next step and the now-position
    x1dot=(x1traj(n+1)-pos2(n,1))/timeinc;
    y1dot=(y1traj(n+1)-pos2(n,2))/timeinc;
    %same with x2 and y2
    x2dot=(x2traj(n+1)-pos4(n,1))/timeinc;
    y2dot=(y2traj(n+1)-pos4(n,2))/timeinc;
    %calculation of forces
    forcel(n, 1)=PROPPL(end,4)*acctraj(n,1)/2+sin(postraj(n,3))*PAR(5,2)/PROPPL(end,3)*acctraj(n,3);
    forcer(n, 1)=PROPPL(end,4)*acctraj(n,1)/2-sin(postraj(n,3))*PAR(5,2)/PROPPL(end,3)*acctraj(n,3);
    forcel(n, 2)=PROPPL(end,4)*acctraj(n,2)/2-cos(postraj(n,3))*PAR(5,2)/PROPPL(end,3)*acctraj(n,3)+PROPPL(end,4)*g/2;
    forcer(n, 2)=PROPPL(end,4)*acctraj(n,2)/2+cos(postraj(n,3))*PAR(5,2)/PROPPL(end,3)*acctraj(n,3)+PROPPL(end,4)*g/2;
    %jacobian of left and right arm
    Jl=[-sin(teta(1))*PROP1(end,3)-sin(teta(1)+teta(2))*PROP2(end,3), -sin(teta(1)+teta(2))*PROP2(end,3); cos(teta(1))*PROP1(end,3)+cos(teta(1)+teta(2))*PROP2(end,3), cos(teta(1)+teta(2))*PROP2(end,3)];
    Jr=[-sin(teta(3))*PROP3(end,3)-sin(teta(3)+teta(4))*PROP4(end,3), -sin(teta(3)+teta(4))*PROP4(end,3); cos(teta(3))*PROP3(end,3)+cos(teta(3)+teta(4))*PROP4(end,3), cos(teta(3)+teta(4))*PROP4(end,3)];
    Jlm=[-sin(teta(1))*PROP1(end,3); cos(teta(1))*PROP1(end,3)];
    Jrm=[-sin(teta(3))*PROP3(end,3); cos(teta(3))*PROP3(end,3)];
    J1=[-sin(teta(1))*PROP1(end,3)*PAR(1,1); cos(teta(1))*PROP1(end,3)*PAR(1,1)];
    J2=[-sin(teta(1))*PROP1(end,3)-sin(teta(1)+teta(2))*PROP2(end,3)*PAR(2,1), -sin(teta(1)+teta(2))*PROP2(end,3)*PAR(2,1); cos(teta(1))*PROP1(end,3)+cos(teta(1)+teta(2))*PROP2(end,3)*PAR(2,1), cos(teta(1)+teta(2))*PROP2(end,3)*PAR(2,1)];
    J3=[-sin(teta(3))*PROP3(end,3)*PAR(3,1); cos(teta(3))*PROP3(end,3)*PAR(3,1)];
    J4=[-sin(teta(3))*PROP3(end,3)-sin(teta(3)+teta(4))*PROP4(end,3)*PAR(4,1), -sin(teta(3)+teta(4))*PROP4(end,3)*PAR(4,1); cos(teta(3))*PROP3(end,3)+cos(teta(3)+teta(4))*PROP4(end,3)*PAR(4,1), cos(teta(3)+teta(4))*PROP4(end,3)*PAR(4,1)];
    %tetadot is the angular velocity at the joints
    tetadotold=tetadot;
    tetadot=[Jl^-1, zeros(2); zeros(2), Jr^-1]*[x1dot;y1dot;x2dot;y2dot];
    tetaddot=(tetadot-tetadotold)/timeinc;
    %Jointspeeds:
    v1old=v1; v2old=v2; v3old=v3; v4old=v4;
    v1=Jlm*tetadot(1);
    v2=Jl*[tetadot(1); tetadot(2)];
    v3=Jrm*tetadot(3);
    v4=Jr*[tetadot(3); tetadot(4)];
    accj1=(v1-v1old)/timeinc;
    accj2=(v2-v2old)/timeinc;
    accj3=(v3-v3old)/timeinc;
    accj4=(v4-v4old)/timeinc;
    accj=[accj1'; accj2'; accj3'; accj4'];
    
    %calculation of the link accelerations
    linkveltraj(n,1:2)=(J1*tetadot(1))';
    linkveltraj(n,3:4)=(J2*[tetadot(1); tetadot(2)])';
    linkveltraj(n,5:6)=(J3*tetadot(3))';
    linkveltraj(n,7:8)=(J4*[tetadot(3); tetadot(4)])';
    if n==1
        accl1=linkveltraj(n, 1:2)'/timeinc;
        accl2=linkveltraj(n, 3:4)'/timeinc;
        accl3=linkveltraj(n, 5:6)'/timeinc;
        accl4=linkveltraj(n, 7:8)'/timeinc;
    else
        accl1=(linkveltraj(n, 1:2)-linkveltraj(n-1, 1:2))'/timeinc;
        accl2=(linkveltraj(n, 3:4)-linkveltraj(n-1, 3:4))'/timeinc;
        accl3=(linkveltraj(n, 5:6)-linkveltraj(n-1, 5:6))'/timeinc;
        accl4=(linkveltraj(n, 7:8)-linkveltraj(n-1, 7:8))'/timeinc;
    end
        
    %calculation of tau
    taul(n,:)=(Jl'*-(forcel(n,:))'+J2'*-PROP2(end,4)*accl2)';
    taur(n,:)=(Jr'*-(forcer(n,:))'+J4'*-PROP4(end,4)*accl4)';
    taul(n,1)=taul(n,1)+(J1'*-PROP1(end,4)*accl1)'-PAR(1,2)*(tetadot(1)-tetadotold(1))/timeinc-PAR(2,2)*(tetadot(1)+tetadot(2)-tetadotold(1)-tetadotold(2))/timeinc;
    taul(n,2)=taul(n,2)-PAR(2,2)*(tetadot(1)+tetadot(2)-tetadotold(1)-tetadotold(2))/timeinc;
    taur(n,1)=taur(n,1)+(J3'*-PROP3(end,4)*accl3)'-PAR(3,2)*(tetadot(3)-tetadotold(3))/timeinc-PAR(4,2)*(tetadot(3)+tetadot(4)-tetadotold(3)-tetadotold(4))/timeinc;
    taur(n,2)=taur(n,2)-PAR(4,2)*(tetadot(3)+tetadot(4)-tetadotold(3)-tetadotold(4))/timeinc;
    
    %the new teta is the old teta plus the joint velocities times the time
    %increment
    
    deform=deformations(accj, forcel(n,:), forcer(n,:), PROP1, PROP2, PROP3, PROP4, PAR, teta, taul(n,:), taur(n,:));
    
    teta=teta+tetadot*timeinc;
end
n=n+1;
 pos2(n,1)=cos(teta(1))*PROP1(end,3)-sin(teta(1))*deform(1,1)+cos(teta(1)+teta(2)+deform(1,2))*PROP2(end,3)-sin(teta(1)+teta(2))*deform(2,1);
 pos2(n,2)=sin(teta(1))*PROP1(end,3)+cos(teta(1))*deform(1,1)+sin(teta(1)+teta(2)+deform(1,2))*PROP2(end,3)+cos(teta(1)+teta(2))*deform(2,1);
 pos4(n,1)=l0+cos(teta(3))*PROP3(end,3)-sin(teta(3))*deform(3,1)+cos(teta(3)+teta(4)+deform(3,2))*PROP4(end,3)-sin(teta(3)+teta(4))*deform(4,1);
 pos4(n,2)=sin(teta(3))*PROP3(end,3)+cos(teta(3))*deform(3,1)+sin(teta(3)+teta(4)+deform(3,2))*PROP4(end,3)+cos(teta(3)+teta(4))*deform(4,1);
    

tetatraj(n,:)=teta';


%calculation finished, here is just output
pos2
pos4
tetatraj
forcel
forcer
taul
taur
plot(pos2(:,1), pos2(:,2))
axis equal
hold on
plot(pos4(:,1), pos4(:,2))
hold off
%plot(pos4(:,1))

end

end

