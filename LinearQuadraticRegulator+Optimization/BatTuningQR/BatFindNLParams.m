% BAT Algorithm to tune Parameters

%Algorithm 
clear all;
close all;
load('MarchData');
% load('B1North_R_T.mat');
% load('B1N_SWT.mat');
% load('B1N_RWT.mat');
% load('B1S_RWT.mat');
% load('B1S_SWT.mat');
% load('B1South_R_T.mat');
% load('B1Amb_TMAX.mat');

global  B1North_R_T B1N_RWT B1N_SWT B1South_R_T B1S_RWT B1S_SWT B1Amb_TMAX    ParameterVals  
global State1 State2 State3  State4

 


pulseFactor = [   1     0     0     0     0     0;     
                  0  1E-5     0     0     0     0;    
                  0     0  1E-2     0     0     0;   
                  0     0     0    10    0     0;  
                  0     0     0     0  0.01     0; 
                  0     0     0     0     0  0.01];

              
freqFactor = [1;
              1E-5;
              1E-2;
              10;
              0.01;
              0.01;
             ];
          
u1= 166.5;
u2= 0.078;
u3= 0.152;
u4= 8635.;
u5= 0.6601;
u6= 0.38;         
%Initialize Bat Algorithm related parameters

PopSize = 20; % no. of bats doing the search
GenSize = 200; % no. of generations of pulses.
Loudness(1:PopSize) = 0.95;
pulseRate0(1:PopSize) = 0.95;
Fmin = zeros(6,1);
Fmax = 60 * freqFactor;
d = 6; % dimensions of the solution - the parameters under scanner
%  ParameterVals = [u1 u2 u3 u4 u5 u6];
% ParameterVals = [99.9000000000000,0.0913237876133860,0.180149637072743,5181,0.396060000000000,0.407461977612594;]
% ParameterVals = [139.860000000000,0.0688203221081286,0.108089782243646,3108.60000000000,0.237636000000000,0.244477186567556;];
% ParameterVals = [195.205103287342,0.0569594248857887,0.0659505524679284,1866.58816485222,0.255448015686383,0.280218237395590;];
% ParameterVals  = [127.080978918424,0.0386950390556416,0.0577084958428220,1886.06051400100,0.236181633197954,0.168130942437354;];
% ParameterVals = [103.682230050956,0.0490007432733878,0.0740819484222587,1434.60488638191,0.304566819805006,0.193074372654089;];
% ParameterVals = [92.7735817934094,0.0539250655491382,0.0906897839401676,953.869671370222,0.409188195522471,0.235660814618166;]
% ParameterVals = [77.2844450452856,0.0728799061493820,0.120133194823972,801.323351949469,0.505840268054139,0.295348998459103;];
ParameterVals  = [60.8862006693394,0.0957077431086720,0.154791328573516,632.668206858105,0.647746646047895,0.177429929211991;];
% Lower limit/bounds/ a vector
Lb= 0.6*ParameterVals;


% Upper limit/bounds/ a vector
Ub= 1.4*ParameterVals;

% Initializing arrays
Q=zeros(PopSize,1);   % Frequency
v=zeros(PopSize,d);   % Velocities

% Initialize the population/solutions
for i=1:PopSize,
  Sol(i,:) = Lb+(Ub-Lb).*rand(1,d);
  Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
  Fitness(i)= GetFitNL(Sol(i,:));
  while ((Fitness(i) < 1) || (Fitness(i)>10))
      Sol(i,:) = Lb+(Ub-Lb).*rand(1,d);
      Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
      Fitness(i)= GetFitNL(Sol(i,:));
  end
i 
Fitness (i)
end

[fmin,I]=max(Fitness);

best=Sol(I,:);
% for i=1:PopSize,
%   if Fitness(i) < 0
%       Sol(i,:) = best;
%   end  
% end
N_iter=0;       % Total number of function evaluations
ischanged(1:PopSize) = 0;
Fprogress = 1;
for t=1:GenSize, 
    t
    
%     if (ischanged == 1)
%         Loudness = Loudness *0.95;
% %         pulseRate0 = pulseRate0 * 1.05;
%         ischanged = 0;
%     end
% Loop over all bats/solutions
        for i=1:PopSize,
          if (ischanged(i) == 1)
              ischanged(i) = 0;
              Loudness (i) = Loudness (i) * 0.9;
              pulseRate0(i) = pulseRate0 (i) * 0.9;
%           else
%               Loudness (i) = Loudness (i) * 1.05;
%               pulseRate0(i) = pulseRate0 (i) * 1.05; 
          end
%               Loudness (i) = Loudness (i) * 0.95;
%               pulseRate0(i) = pulseRate0 (i) * 0.9;    

          for countx = 1:d
            Q(i,countx)=Fmin(countx)+(Fmin(countx)-Fmax(countx))*(rand);
            v(i,countx)=v(i,countx)+(Sol(i,countx)-best(countx))*Q(i,countx);
          end
%           v(i,:)=v(i,:)+(Sol(i,:)-best)*Q(i,:);
          S(i,:)=Sol(i,:)+v(i,:);
          % Apply simple bounds/limits
          Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
          % Pulse rate
          if rand>pulseRate0(i)
          % The factor 0.001 limits the step sizes of random walks 
              S(i,:)=best+ (0.1*pulseFactor * randn(1,d)')';
          end

     % Evaluate new solutions
           S(i,:) = simplebounds(S(i,:),Lb,Ub);
           Fnew= GetFitNL(S(i,:));
     % Update if the solution improves, or not too loud
           if (Fnew >= Fitness(i)) & (rand<Loudness (i)) ,
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
           end

          % Update the current best solution
          if Fnew>fmin,
                best=S(i,:);
                fmin=Fnew;
                FitnessStore (Fprogress) = Fnew;
                BestfitParam(Fprogress,1:6) = best;
                Fprogress = Fprogress + 1;
                ischanged(i) = 1;
                best
                ParameterVals
        
                Fnew
                if (Fnew > 70)
                figure;
                subplot(2,2,1);
                plot( B1N_RWT,'k');
                hold on;
                plot(State1);
                subplot(2,2,2);
                plot( B1S_RWT,'k');
                hold on;
                plot(State2);
                subplot(2,2,3);
                plot( B1North_R_T,'k');
                hold on;
                plot(State3);
                subplot(2,2,4);
                plot( B1South_R_T,'k');
                hold on;
                plot(State4);
                
                
             
                end
          end
        end
        
%         for i=1:PopSize,
%             Sol(i,:) = best;
%         end
        N_iter=N_iter+PopSize;
end


% save('BatFitPatientNL1','BestfitParam','FitnessStore');
