function Performance = LQR_QR_Tuning(xinp)
global theta_swn theta_sws theta_o

load('TempDetail');
Rval = zeros(3,3);

Rval(1,1) = xinp(1);
Rval(2,2) = xinp(2);
Rval(3,3) = xinp(3);

Rval1 = zeros(3,3);
Rval1(1,1) = xinp(4);
Rval1(2,2) = xinp(5);
Rval1(3,3) = xinp(6);

Qval = zeros(4,6);

Qval(1,3) = xinp(7);
Qval(2,4) = xinp(8);
Qval(3,5) = xinp(9);
Qval(4,6) = xinp(10);


u1=50.3172;
u2=0.1199;
u3=0.1920;
u4=515.0689;
u5=0.7980;
u6=0.1315;


A=[ -( (1/(u1*u2)) + (1/(u1*u3)) ) , 0 , (1/(u1*u2)) , 0 ;
    0 , -( (1/(u1*u2)) + (1/(u1*u3)) ), 0 ,(1/(u1*u2)) ;
 (1/(u4*u2)) , 0  ,-( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) ) ,(1/(u4*u6));
    0 ,(1/(u4*u2)) ,(1/(u4*u6)),-( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) )];


B=[(1/(u1*u3)),0,0;
   0,(1/(u1*u3)),0;
   0,0,(1/(u4*u5));
   0,0,(1/(u4*u5))];
% 
C=[0,0,1,0;
    0,0,0,1];

D=0;

sys=ss(A,B,C,D);%creates a state-space model object representing the continuous-time state-space model

sysd = c2d(sys,1);%Convert model from continuous to discrete time
A = sysd.a;
B = sysd.b;
C = sysd.c;
D = sysd.d;

R = Rval;

R1 = Rval1;

[i,j] = size(A);
% Q = 1*eye(3);
% Q = 1*eye(4);
Q = 1*eye(4);
% Q(1,1) = 0.0001*Q(1,1);
% Q(2,2) = 0.0001*Q(2,2);


%acompZero = 0;
acompZero = zeros(1,2);
for count = 1:(j-1)
%     acompZero = vertcat(acompZero,0);
acompZero = vertcat(acompZero,zeros(1,2));
end


Anew = horzcat(A,acompZero);
[i,j] = size(C*A);

for count = 1:(i-1)
    acompZero = horzcat(acompZero,acompZero);
end

if (i == 1)
    icomp = 1;
elseif (i == 2)
    icomp = [1 0;0 1];
elseif (i == 3)
    icomp = [1 0 0;0 1 0;0 0 1];
end

Anewrow2 = [-(C*A) icomp];
Anew =  vertcat(Anew,Anewrow2);
Bnew = vertcat(B, -(C*B));

Qcomp = Qval;



Qnew = Qcomp'*Q*Qcomp;

[Qnewt,L,G] = dare(Anew,Bnew,Qnew,R);

[i,j] = size(Qnew);

for count = 1:i
    for count1 = 1:j
        QnewComp(count,count1) = 0;
    end
end 

predictionHorizon = 6;

dump=zeros(6,6);

GammaX=vertcat(horzcat(Qnew, dump, dump, dump, dump,dump), ...       
               horzcat(dump, Qnew, dump, dump, dump,dump), ...   
               horzcat(dump, dump, Qnew, dump, dump,dump), ...
               horzcat(dump, dump, dump, Qnew, dump,dump), ...
               horzcat(dump, dump, dump, dump, Qnew,dump), ...
               horzcat(dump, dump, dump, dump, dump,Qnewt));
           
dump1 = zeros(3,3); 

GammaU=vertcat(horzcat(R , dump1, dump1, dump1, dump1, dump1),...
               horzcat(dump1,  R, dump1, dump1, dump1, dump1),...
               horzcat(dump1, dump1,  R, dump1, dump1, dump1),...
               horzcat(dump1, dump1, dump1, R,  dump1, dump1),...
               horzcat(dump1, dump1, dump1, dump1,  R, dump1),...
               horzcat(dump1, dump1, dump1, dump1,dump1,   R));
    
           
           
dumpx = zeros(3,3); 

GammaU1=vertcat(horzcat(R1 , dump1, dumpx, dumpx, dumpx, dumpx),...
               horzcat(dumpx,  R1, dumpx, dumpx, dumpx, dumpx),...
               horzcat(dumpx, dumpx,  R1, dumpx, dumpx, dumpx),...
               horzcat(dumpx, dumpx, dumpx, R1,  dumpx, dumpx),...
               horzcat(dumpx, dumpx, dumpx, dumpx,  R1, dumpx),...
               horzcat(dumpx, dumpx, dumpx, dumpx,dumpx,   R1));          
           
dump2 = zeros(3,6);
GammacP=[1 0 0 0 0 0;
           0 1 0 0 0 0
           0 0 0 0 0 0];
      
   GammaC=vertcat(horzcat(GammacP, dump2, dump2, dump2, dump2, dump2), ...
                  horzcat(dump2, GammacP, dump2, dump2, dump2, dump2), ...       
                  horzcat(dump2, dump2, GammacP, dump2, dump2, dump2), ...   
                  horzcat(dump2, dump2, dump2, GammacP, dump2, dump2), ...
                  horzcat(dump2, dump2, dump2, dump2, GammacP, dump2), ...
                   horzcat(dump2, dump2, dump2, dump2, dump2, GammacP));
%                horzcat(dump2, dump2, dump2, dump2, dump2, dump2,GammacP));

[i,j] = size(Anew);

for count = 1:i
    for count1 = 1:j
        if (count == count1)
            I(count,count1) = 1;
        else
            I(count,count1) = 0;
        end
    end
end

for count = 1:(predictionHorizon )%Changed
    if(count == 1)
        Sx = Anew;
    else
        Sx = vertcat(Sx,Anew^(count-1));
    end
end
[i,j] = size(Bnew);

for count = 1:i
    for count1 = 1:j
        BnewComp(count,count1) = 0;
    end
end

for count = 1:predictionHorizon
    if(count == 1)
        Su = Bnew;
        for count1 = 2:predictionHorizon
            Su = horzcat(Su,BnewComp);
        end

    else
        currCount = count - 2;
        SuTemp = 0;
        for  count1= 1:predictionHorizon
            if (currCount == (count-2))
            SuTemp = Anew^(currCount)*Bnew;
            currCount = currCount -1;
            elseif (currCount >= 0)
            
                SuTemp = horzcat(SuTemp,Anew^(currCount)*Bnew);
            
            currCount = currCount -1;
            else
             SuTemp = horzcat(SuTemp,BnewComp);
          
            end
        end
        
        Su = vertcat(Su,SuTemp);
        
    end
end

H = (GammaU + Su'*GammaC'*GammaU1*GammaC*Su - GammaU1*GammaC*Su - Su'*GammaC'*GammaU1)+ (Su'*GammaX* Su);
g = (Sx'*GammaC'*GammaU1*GammaC*Su - Sx'*GammaC'*GammaU1 + Sx'*GammaX*Su)';


ParameterVals  = [50.3172077776522,0.119911737440304,0.192017197106273,515.068934326141,0.798016930154889,0.131480458729525;];

k = 1;

X(1,1:6)=0;
x(1,1:6) = 0;
x_ss(1:4) = [B1N_RWT(1,1); B1S_RWT(1,1);20;20];
x_ss(5:10) = ParameterVals(1:6);

theta_swn = B1N_SWT(1,1);               
theta_sws = B1S_SWT(1,1);
 
theta_o = B1Amb_TMAX (1,1);
% [t,xx] = ode45('Tuning_ODE',[0 1000],x_ss);
% [i,j] = size(t);
% x_ss(1:4) = xx(2000,1:4);
% x_ss(5:10) = ParameterVals(1:6);

F=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
   0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
   0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];                  
ExptLength = 2500;

State1 = zeros(ExptLength,1);
State2 = zeros(ExptLength,1);
State3 = zeros(ExptLength,1);
State4 = zeros(ExptLength,1);

[i,j]=size(B1Amb_TMAX);

diff_amb(1)=0;
for count=1:(j-1)
    diff_amb(count+1)=B1Amb_TMAX(count+1)-B1Amb_TMAX(count);
end

        

 ref1 = 20;
 ref2 = 20;

countx = 1;
 diffambTemp = 0;

while(k<=ExptLength)

Ux(k,1:3)= -F*inv(H)*g*x(k,1:6)';


if (theta_swn + Ux(k,1) > 70)
    theta_swn = 70;
    Ux(k,1)= 0;
elseif (theta_swn + Ux(k,1) <15)
    theta_swn = 15;
    Ux(k,1)= 0;
else
    theta_swn = theta_swn + Ux(k,1) ;
end

if (theta_sws + Ux(k,2) > 70)
    theta_sws = 70;
    Ux(k,2) = 0;
elseif (theta_sws + Ux(k,2) <15)
    theta_sws = 15;
    Ux(k,2) = 0;
else
    theta_sws = theta_sws + Ux(k,2) ;
end

swn(k) = theta_swn;
sws(k) = theta_sws;
amb(k) = theta_o;


if(k==2000)
    ref1=25;
    ref2=25;

end

if (mod(k,18) == 0)
%  theta_o = B1Amb_TMAX(1,countx);
  theta_o = B1Amb_TMAXJan(1,countx);
%  theta_o = B1Amb_TMAXMay(1,countx);
 countx = countx + 1;
 if (countx >1)
%  diffambTemp = B1Amb_TMAXMay(1,countx) - B1Amb_TMAXMay(1,countx-1);
diffambTemp = B1Amb_TMAXJan(1,countx) - B1Amb_TMAXJan(1,countx-1);
% diffambTemp = B1Amb_TMAX(1,countx) - B1Amb_TMAX(1,countx-1);
 else
 diffambTemp = 0;
 end
end
Ux(k,3)= (1/18) * diffambTemp;

[t,xx] = ode45('Tuning_ODE',[0 1],x_ss);
[i,j] = size(t);
x_ss(1:4) = xx(i,1:4);
x_ss(5:10) = ParameterVals(1:6);


x(k+1,1:4)= A*x(k,1:4)' + B*Ux(k,1:3)';
State1(k) =  x_ss(1);
State2(k) =  x_ss(2);
State3(k) =  x_ss(3);
State4(k) =  x_ss(4);


x(k+1,5) = ref1 - x_ss(3);
x(k+1,6) = ref2 - x_ss(4);
ref1Data(k)  = ref1;
amb(count)= theta_o;

 k = k + 1;
end

% Emmeasure = (sum((swn(2000:2200)' - State1(2000:2200))) + sum((sws(2000:2200)' - State2(2000:2200))))/18';
% Emmeasure =   (sum((swn(2000:2250)' - State1(2000:2250))) + sum((sws(2000:2250)' - State2(2000:2250))));
Emmeasure =(sum((swn(2000:2300)' - State1(2000:2300))) + sum((sws(2000:2300)' - State2(2000:2300))));
msemeasure = mse(ref1Data(1000:ExptLength)' - State3(1000:ExptLength));

Performance = [1/Emmeasure msemeasure];