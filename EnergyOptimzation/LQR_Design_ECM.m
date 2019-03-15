
%Discretizing system state A , B , C , D
clc
clear all
load('TempDetail.mat')
load('Noise.mat')
% parameters

u1=50.3172;
u2=0.1199;
u3=0.1920;
u4=515.0689;
u5=0.7980;
u6=0.1315;

% a11= -( (1/(u1*u2)) + (1/(u1*u3)) ) ;
% a13=  (1/(u1*u2));
% a22= -( (1/(u1*u2)) + (1/(u1*u3)) ) ;
% a24=  (1/(u1*u2));
% a31= (1/(u4*u2));
% a33= -( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) );
% a34= (1/(u4*u6));
% a42= (1/(u4*u2));
% a43= (1/(u4*u6));
% a44= -( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) );



A=[ -( (1/(u1*u2)) + (1/(u1*u3)) ) , 0 , (1/(u1*u2)) , 0 ;
    0 , -( (1/(u1*u2)) + (1/(u1*u3)) ), 0 ,(1/(u1*u2)) ;
 (1/(u4*u2)) , 0  ,-( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) ) ,(1/(u4*u6));
    0 ,(1/(u4*u2)) ,(1/(u4*u6)),-( (1/(u4*u5)) + (1/(u4*u6)) + (1/(u4*u2)) )];


% B=[(1/(u1*u3)),0;
%    0,(1/(u1*u3));
%    0,0;
%    0,0];


B=[(1/(u1*u3)),0,0;
   0,(1/(u1*u3)),0;
   0,0,(1/(u4*u5));
   0,0,(1/(u4*u5))];
% 
% C=[0,0,1,0;
%     0,0,0,1];

C=[ 1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];

% K = [0;
%      0;
%      (1/(u4*u5));
%      (1/(u4*u5));];

D=0;

sys=ss(A,B,C,D);%creates a state-space model object representing the continuous-time state-space model

sysd = c2d(sys,1);%Convert model from continuous to discrete time
A = sysd.a;
B = sysd.b;
C = sysd.c;
D = sysd.d;

sys=ss(A,B,C,D);

sysd=c2d(sys,1); % Discrete State Space model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(Section  for devloping Anew, Bnew matrix%%%%%%%%%%%%%%%%%%%%%%

% sys=ss(A,B,C,D); %%% continuous time state space design
% 
% sysd=c2d(sys,1);%%% Discretization of  state space

% Updating discretizing State Space

 A=sysd.a;
 B=sysd.b;
 C=sysd.c;
 D=sysd.d;
 
 %%% Creating new A matrix pertaining to set point tracking
 
 [row_A,column_A]=size(A);
 
 [row_CA,column_CA]=size(C*A); 
 
 Temp=zeros(row_A,row_CA); %%% creating zero padding to horzcat in A
 
 Anew=horzcat(A,Temp);
 
 Anewrow=horzcat(-(C*A),eye(row_CA));
 
 Anew=vertcat(Anew,Anewrow); % creating final Anew square matrix
 
 %%% Creating new B matrix pertaining to set point tracking
 
 Bnew=vertcat(B,-(C*B));
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%Change coefficients of Q and R matrix aptly to obtain desired response%%%%%%
 
 [row_Anew,column_Anew]=size(Anew); % determining size of A to create Q matrix
 
 Q=1*eye(row_Anew); % Definig Q matrix which must match to the size of Anew and Tuning parameter for minimizing error in states
 
 [row_B,column_B]=size(B); % Determining  no of input from cloumn size to define R 
 
 R=0.5*eye(column_B); % Tuning parameter for input signal
 
 Qt=dare(Anew,Bnew,Q,R); % algebric riccati equation solved to get new matrix Qt to take end diagonal position in gammax
 
 
 % Define  prediction horizon
 
 Prediction=6;
 

% Defining New state space matrix based on prediction horizon 

% New A matrix named Sx used for  designing Hessian and gradient

I=eye(row_Anew); % size of identity matrix must match with suared matrix Anew

Sx=I;

    for i=1:Prediction;
        Sx=vertcat(Sx,Anew^i);
    end
   
% New B matrix named Su used for designing Hessian and gradient

[row_Bnew,column_Bnew]= size(Bnew);

    for i=1:Prediction;
        Su1=zeros(i*row_Anew,column_Bnew);
        for j=i:Prediction;
            Su1=vertcat(Su1,Anew^(j-i)*Bnew);
        end
        if (i==1)
            Su=Su1;
        else
            Su=horzcat(Su,Su1);
        end
    end
 
% New Tuning matrix Q named Gammax for designing input u
    
    GammaX=Q;
    
    for i=1:Prediction;
        if(i==Prediction)
            GammaX=blkdiag(GammaX,Qt);
        else
            GammaX=blkdiag(GammaX,Q);
        end
    end
    
 % New Tuing matrix R named Gammau for designing u  
    
    GammaU=R;
    
    for i=1:(Prediction-1);
        GammaU=blkdiag(GammaU,R);
    end
    
           
           
dump2 = zeros(2,6);
GammacP=2*[1 0 0 0 0 0;
       0 1 0 0 0 0];
      
   GammaC=vertcat(horzcat(GammacP, dump2, dump2, dump2, dump2, dump2,dump2), ...
               horzcat(dump2, GammacP, dump2, dump2, dump2, dump2,dump2), ...       
               horzcat(dump2, dump2, GammacP, dump2, dump2, dump2,dump2), ...   
               horzcat(dump2, dump2, dump2, GammacP, dump2, dump2,dump2), ...
               horzcat(dump2, dump2, dump2, dump2, GammacP, dump2,dump2), ...
               horzcat(dump2, dump2, dump2, dump2, dump2, GammacP,dump2));
%                horzcat(dump2, dump2, dump2, dump2, dump2, dump2,GammacP));


% H=-(Sx'*GammaC'*GammaU)-(Sx'*GammaC'*GammaU*GammaC*Su)+(Su'*GammaX*Su)+(GammaU)+(GammaU*GammaC*Su)...
%     +(Su'*GammaC'*GammaU)+(Su'*GammaC'*GammaU*GammaC*Su);
%  
% g=(Sx'*GammaX*Sx)+(Su'*GammaX*Sx)-(GammaU*GammaC*Sx)+(Sx'*GammaC'*GammaU*GammaC*Sx)-(Su'*GammaC'*GammaU*GammaC*Sx);

H=(Su'*GammaX*Su)+(GammaU)+(GammaU*GammaC*Su)+(Su'*GammaC'*GammaU)+(Su'*GammaC'*GammaU*GammaC*Su);
g=(Su'*GammaX*Sx)-(GammaU*GammaC*Sx)-(Su'*GammaC'*GammaU*GammaC*Sx);


         
% H = Su'*GammaX*Su + GammaU;
% g = Su'*GammaX*Sx;

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
global  theta_swn theta_sws theta_o 
% ExptLength = 300;
ParameterVals  = [50.3172077776522,0.119911737440304,0.192017197106273,515.068934326141,0.798016930154889,0.131480458729525;];

k = 1;

iterations=3;
X(1,1:6)=0;
x(1,1:6) = 0;
x_ss(1:4) = [B1N_RWT(1,1); B1S_RWT(1,1);B1North_R_T(1,1);B1South_R_T(1,1)];
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
F=[1 0 0 0 0 0 0 0 0 0 0 0 ; ...
                      0 1 0 0 0 0 0 0 0 0 0 0 ];                  
ExptLength = 5500;

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

countx = 1;
 diffambTemp = 0;
while(k<=ExptLength)

Ux(k,1:2)= -F*inv(H)*g*x(k,1:6)';

if (theta_swn + Ux(k,1) > 70)
    theta_swn = 70;
elseif (theta_swn + Ux(k,1) < 15)
    theta_swn = 15;
else
    theta_swn = theta_swn + Ux(k,1);
end
  
if (theta_sws + Ux(k,2) > 70)
    theta_sws = 70;
elseif (theta_sws + Ux(k,2) < 15)
    theta_sws = 15;
else
    theta_sws = theta_sws + Ux(k,2);
end


           
% theta_swn = theta_swn + Ux(k,1) ;
% theta_sws = theta_sws + Ux(k,2) ;
% theta_o = theta_o + Ux(k,3);
swn(k) = theta_swn;
sws(k) = theta_sws;
amb(k) = theta_o;

% if(k==1000)
%     ref1=25;
%     ref2=25;
% elseif(k==2000)
%     ref1=22;
%     ref2=22;
% elseif(k==3000)
%     ref1=26;
%     ref2=26;
% elseif(k==4000)
%     ref1=23;
%     ref2=23;
% end

if (k==2000)
    ref1=25;
    ref2=25;
end

if (mod(k,18) == 0)
 theta_o = B1Amb_TMAX(1,countx);
 countx = countx + 1;
 if (countx >1)
%  diffambTemp = B1Amb_TMAXMay(1,countx) - B1Amb_TMAXMay(1,countx-1);
% diffambTemp = B1Amb_TMAXJan(1,countx) - B1Amb_TMAXJan(1,countx-1);
diffambTemp = B1Amb_TMAX(1,countx) - B1Amb_TMAX(1,countx-1);
 else
 diffambTemp = 0;
 end
end
% Ux(k,3)= (1/18) * diffambTemp;
Ux(k,3)= 0;
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

x(k+1,1:4)= A*x(k,1:4)' + B*Ux(k,1:2)';
State1(k) =  x_ss(1);
State2(k) =  x_ss(2);
State3(k) =  x_ss(3);
State4(k) =  x_ss(4);

x(k+1,5) = ref1 - x_ss(3)+Feed_Noise(k,1);
x(k+1,6) = ref2 - x_ss(4)+Feed_Noise(k,1);
ref1Data(k)  = ref1;
amb(count)= theta_o;

 k = k + 1;
end

figure;
plot(State3, 'r');
hold on;
% plot(NBT(:,2));
% Emmeasure = sum(swn - State1) + sum(sws - State2);
%  Emmeasure = sum(abs(swn' - State1)) + sum(abs(sws' - State2));
