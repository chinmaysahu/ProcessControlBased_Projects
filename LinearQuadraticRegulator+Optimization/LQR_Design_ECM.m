
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
C=[0,0,1,0;
    0,0,0,1];

% C=[ 1,0,0,0;
%     0,1,0,0;
%     0,0,1,0;
%     0,0,0,1];

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
% K = [0.0196;0.0196;0.0395;0.0395];
% 
% sysAmb = ss(A,BAmb,C,0);
% sysdAmb = c2d(sysAmb,1);%Convert model from continuous to discrete time
% Aamb = sysdAmb.a;
% Bamb = sysdAmb.b;
% Camb = sysdAmb.c;
% Damb = sysdAmb.d;

% R = 1;
%   R=0.0001*[1 0 0;
%      0 1 0;
%      0 0 1];

R=0.1*[1 0 0;
       0 1 0
       0 0 1];    
   
R1=0.5*[1 0 0;
       0 1 0
       0 0 1]; 

[i,j] = size(A);
% Q = 1*eye(3);
% Q = 1*eye(4);
Q = 0.01*eye(4);
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


[i,j] = size(Q);

if (i == 1)
    Qcomp = [0 1];
elseif (i == 2)
    Qcomp = [0 1 0;
             0 0 1];
elseif (i == 3)
    Qcomp = [0 1 0 0;
             0 0 1 0;
             0 0 0 1];
elseif (i == 4)
    Qcomp = [0 0 1 0 0 0 ;
             0 0 0 1 0 0 ;
             0 0 0 0 20 0 ;
             0 0 0 0 0 20 ;
            ];
  elseif (i == 5)
    Qcomp = [0 0 1 0 0 0 ;
             0 0 0 1 0 0 ;
             0 0 0 0 20 0 ;
             0 0 0 0 0 20 ;
            ];       
            

end



Qnew = Qcomp'*Q*Qcomp;

[Qnewt,L,G] = dare(Anew,Bnew,Qnew,R);

[i,j] = size(Qnew);

for count = 1:i
    for count1 = 1:j
        QnewComp(count,count1) = 0;
    end
end 

predictionHorizon = 6;

% for count4 = 1:(predictionHorizon + 1) 
%     for count5 = 1:(predictionHorizon + 1)
%         if ((count4 == count5))
%             if (count5 < (predictionHorizon + 1))
%                 UpdMatrix = Qnew;
%                 Gammau(count4,count5) = R;
%             elseif(count5 == (predictionHorizon+1))
%                 UpdMatrix = Qnewt;
%                
%             end
%            
%         else
%             if ((count5 < (predictionHorizon + 1)) && (count4 < (predictionHorizon + 1)))
%              Gammau(count4,count5) = 0;
%             end
%             UpdMatrix = QnewComp;
%             
%         end
%         
%         for count = 1:i
%             for count1 = 1:j
%                 Gammax((count4 -1)* i+count,(count5 -1)* j+ count1) = UpdMatrix(count,count1);
%             end
%         end
%     end
% end

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

countx = 1;
 diffambTemp = 0;
while(k<=ExptLength)

Ux(k,1:3)= -F*inv(H)*g*x(k,1:6)';


if (theta_swn + Ux(k,1) > 60)
    theta_swn = 60;
    Ux(k,1)= 0;
elseif (theta_swn + Ux(k,1) <15)
    theta_swn = 15;
    Ux(k,1)= 0;
else
    theta_swn = theta_swn + Ux(k,1) ;
end

if (theta_sws + Ux(k,2) > 60)
    theta_sws = 60;
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

if(k==2000)
    ref1=25;
    ref2=25;
% elseif(k==2000)
%     ref1=22;
%     ref2=22;
% elseif(k==3000)
%     ref1=26;
%     ref2=26;
% elseif(k==4000)
%     ref1=23;
%     ref2=23;
end


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
 theta_o = B1Amb_TMAXMay(1,countx);
%   theta_o = B1Amb_TMAXJan(1,countx);
%  theta_o = B1Amb_TMAXMay(1,countx);
 countx = countx + 1;
 if (countx >1)
 diffambTemp = B1Amb_TMAXMay(1,countx) - B1Amb_TMAXMay(1,countx-1);
% diffambTemp = B1Amb_TMAXJan(1,countx) - B1Amb_TMAXJan(1,countx-1);
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

x(k+1,5) = ref1 - x_ss(3);
x(k+1,6) = ref2 - x_ss(4);
ref1Data(k)  = ref1;
amb(count)= theta_o;

 k = k + 1;
end

figure;
plot(State3, 'r');
hold on;
% plot(NBT(:,2));
% Emmeasure = sum(swn - State1) + sum(sws - State2);
 Emmeasure = sum(abs(swn(1000:18:5000)' - State1(1000:18:5000))) + sum(abs(sws(1000:18:5000)' - State2(1000:18:5000)))
% Emmeasure = sum(swn' - State1) + sum(sws' - State2)