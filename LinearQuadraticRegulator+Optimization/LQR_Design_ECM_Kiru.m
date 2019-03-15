
%Discretizing system state A , B , C , D
clc
clear all
close all
% parameters

u1=50.3172;
u2=0.1199;
u3=0.1920;
u4=515.0689;
u5=0.7980;
u6=0.1315;

% Params_Bat = [6.23838337751987,4.08231143093526,2.65528420692094,3.69013962891634,3.54927049982546,3.04233374621455,1.30789351389086,1.39506660757733,14.9255287416866,9.04766383710071;];
% Params_Bat = [8.41748435115720,5.11059675876186,0.384694552176183,2.33877330020370,0.354969777919614,3.01047791904419,1.32375374137403,2.39160040374453,27.6198898953542,10.3285700566350;];
Params_Bat(1:10) = 1; 
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

[i,j] = size(A);


Q = 1*eye(4);
Qcomp = zeros(4,6);
Qcomp(1,3) = Params_Bat(7);
Qcomp(2,4) = Params_Bat(8);
Qcomp(3,5) = Params_Bat(9);
Qcomp(4,6) = Params_Bat(10);

R = 5* eye(3);
R(1,1) = Params_Bat(1);
R(2,2) = Params_Bat(2);
R(3,3) = Params_Bat(3);

R1 = 3* eye(3);
R1(1,1) = Params_Bat(4);
R1(2,2) = Params_Bat(5);
R1(3,3) = Params_Bat(6);

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


% [i,j] = size(Q);
% 
% if (i == 1)
%     Qcomp = [0 1];
% elseif (i == 2)
%     Qcomp = [0 1 0;
%              0 0 1];
% elseif (i == 3)
%     Qcomp = [0 1 0 0;
%              0 0 1 0;
%              0 0 0 1];
% elseif (i == 4)
%     Qcomp = [0 0 1 0 0 0 ;
%              0 0 0 1 0 0 ;
%              0 0 0 0 10 0 ;
%              0 0 0 0 0 10 ;
%             ];
%   elseif (i == 5)
%     Qcomp = [0 0 1 0 0 0 ;
%              0 0 0 1 0 0 ;
%              0 0 0 0 15 0 ;
%              0 0 0 0 0 15 ;
%             ];       
%             
% 
% end



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

 


%closed loop control

load('B1North_R_T.mat');
load('B1N_SWT.mat');
load('B1N_RWT.mat');
load('B1S_RWT.mat');
load('B1S_SWT.mat');
load('B1South_R_T.mat');
load('B1Amb_TMAX.mat');
load('B1AmbTMAX_MAY.mat');
load('B1AmbTMAXJan.mat');
 load('Feed_Noise.mat')
global  theta_swn theta_sws theta_o 
% ExptLength = 300;
ParameterVals  = [50.3172077776522,0.119911737440304,0.192017197106273,515.068934326141,0.798016930154889,0.131480458729525;];

k = 1;

iterations=3;
X(1,1:6)=0;
x(1,1:6) = 0;
% x_ss(1:4) = [B1N_RWT(1,1); B1S_RWT(1,1);B1North_R_T(1,1);B1South_R_T(1,1)];
x_ss(1:4) = [B1N_RWT(1,1); B1S_RWT(1,1);20;20];
x_ss(5:10) = ParameterVals(1:6);
x_ss(11) = 0;
theta_swn = B1N_SWT(1,1);               
theta_sws = B1S_SWT(1,1);
 
theta_o = B1Amb_TMAX (1,1);
[t,xx] = ode45('Tuning_ODE',[0 1000],x_ss);
[i,j] = size(t);
x_ss(1:4) = xx(2000,1:4);
x_ss(5:10) = ParameterVals(1:6);
x_ss(11) = 0;
% F=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
%                       0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;...
%                       0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
F=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
   0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
   0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];                  
ExptLength = 5500;
% ExptLength = 320;
State1 = zeros(ExptLength,1);
State2 = zeros(ExptLength,1);
State3 = zeros(ExptLength,1);
State4 = zeros(ExptLength,1);

[i,j]=size(B1Amb_TMAX);

diff_amb(1)=0;
for count=1:(j-1)
    diff_amb(count+1)=B1Amb_TMAX(count+1)-B1Amb_TMAX(count);
%       diff_amb (count+1) = 0;
end


                  
% ref1 = x_ss(3);
% ref2 = x_ss(4);
 
 ref1 = 20;
 ref2 = 20;
 RefStd = 20;
countx = 1;
 diffambTemp = 0;
 ChangeRef = 0;
 ChangeRefCount = 1;
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

% theta_sws = theta_sws + Ux(k,2) ;
% theta_o = theta_o + Ux(k,3);

    

swn(k) = theta_swn;
sws(k) = theta_sws;
amb(k) = theta_o;

if (k == 1)

elseif(k==2000)
    ref1=25;
    ref2=25;
    RefStd = 25;

end

% if (ChangeRef == 1)
%     if (ref1 ~= RefStd)
%     ref1 = RefStd;
%     ref2 = RefStd;
%     ChangeRef = 0;
%     ChangeRefCount = 0;
%     else
%         if (ref1 > theta_o)
%         ref1 = ref1*0.95;
%         ref2 = ref1;
%         elseif (ref1 < theta_o)
%         ref1 = ref1 * 1.05;
%         ref2 = ref1;
%         end 
%     ChangeRef = 0;
%     ChangeRefCount = 0;
%     end
% else
%     if ChangeRefCount < 60
%         ChangeRefCount = ChangeRefCount + 1;
%     else
%         ChangeRef = 1;
%     end
% end
%         
    
% 
% if(k==100)
%     ref1=25;
%     ref2=25;
% % elseif(k==150)
% %     ref1=22;
% %     ref2=22;
% elseif(k==200)
%     ref1=23;
%     ref2=23;
% % elseif(k==250)
% %     ref1=23;
% %     ref2=23;
% end


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
% Ux(k,3)=  diffambTemp;
%  theta_o = B1Amb_TMAXJan(1,countx);
% Ux(k,3)= 0;
% % elseif (k > 300)
% % theta_o = B1Amb_TMAX (1,k - 300);
% % Ux(k,3)=diff_amb(k - 300);
% % else
% % theta_o = B1Amb_TMAX (1,k);
% % Ux(k,3)=diff_amb(k);
% % Ux(k,3) = 0;
% % end
%  
% if (k == 350)
% ref1 = 25;
% ref2 = 25; 
% end

% theta_o = 25;

[t,xx] = ode45('Tuning_ODE',[0 1],x_ss);
[i,j] = size(t);
x_ss(1:4) = xx(i,1:4);
x_ss(5:10) = ParameterVals(1:6);
x_ss(11) = 0;

x(k+1,1:4)= A*x(k,1:4)' + B*Ux(k,1:3)';
State1(k) =  x_ss(1);
State2(k) =  x_ss(2);
State3(k) =  x_ss(3);
State4(k) =  x_ss(4);


x(k+1,5) = ref1 - x_ss(3)+0.5*Feed_Noise(k,1);
x(k+1,6) = ref2 - x_ss(4)+0.5*Feed_Noise(k,1);
ref1Data(k)  = ref1;
amb(count)= theta_o;
 mpMPC_SH_SWN(k) = theta_swn;
 mpMPC_SH_SWS(k) = theta_sws;
 k = k + 1;
end

figure;
plot(State3, 'r');
hold on;
 Emmeasure = (sum((swn(2000:2200)' - State1(2000:2200))) + sum((sws(2000:2200)' - State2(2000:2200))))/18
 
 HDD = sum(ref1Data(1000:5000) - amb(1000:5000))/18
 mse(ref1Data(1000:5000)' - State3(1000:5000))
 ise = 0;
 for count = 1000:5000
ise = ise + (abs(ref1Data(count) - State3(count)))^2;
 end
ise

iae = sum(abs(ref1Data(1000:5000)' - State3(1000:5000)))

SH_RT_N = State3;
SH_RT_S = State4;
SH_RWT_N = State1;
SH_RWT_S = State2;


% save('mpMPCSH','SH_RT_N','SH_RT_S','SH_RWT_N','SH_RWT_S','amb','ref1Data','mpMPC_SH_SWN','mpMPC_SH_SWS');