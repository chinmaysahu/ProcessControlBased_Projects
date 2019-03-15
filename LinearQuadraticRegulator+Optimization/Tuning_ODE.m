
function xdot= Tuning_ODE(t,x)

%Initialize Global values here

global theta_swn theta_sws theta_o ambcompensate




x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
u1 = x(5);
u2 = x(6);
u3 = x(7);
u4 = x(8);
u5 = x(9);
u6 = x(10);



%ODE 

x1dot= (1/(u1*u2))*(x3-x1) + (1/(u1*u3))*(theta_swn-x1); %Temp.of return water north.
x2dot= (1/(u1*u2))*(x4-x2) + (1/(u1*u3))*(theta_sws-x2); %Temp.of return water south.
x3dot= (1/(u4*u5))*(theta_o-x3) +  (1/(u4*u6))*(x4-x3) + (1/(u4*u2))*(x1-x3); %Temp. of north block
x4dot= (1/(u4*u5))*(theta_o-x4) + (1/(u4*u6))*(x3-x4) + (1/(u4*u2))*(x2-x4); %Temp. of south block
x5dot = (1/(u4*u5))*(theta_o-x3);

xdot=[x1dot;x2dot;x3dot;x4dot;u1;u2;u3;u4;u5;u6;x5dot];

