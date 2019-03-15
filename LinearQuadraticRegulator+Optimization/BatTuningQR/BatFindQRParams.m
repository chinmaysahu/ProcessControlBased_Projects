% BAT Algorithm to tune Parameters

%Algorithm 
clear all;
close all;

pulseFactor = zeros(10,10);
pulseFactor(1,1) = 0.01;
pulseFactor(2,2) = 0.01;
pulseFactor(3,3) = 0.01;
pulseFactor(4,4) = 0.01;
pulseFactor(5,5) = 0.01;
pulseFactor(6,6) = 0.01;
pulseFactor(7,7) = 0.01;
pulseFactor(8,8) = 0.01;
pulseFactor(9,9) = 0.01;
pulseFactor(10,10) = 0.01;


freqFactor(1) = 0.01;
freqFactor(2) = 0.01;
freqFactor(3) = 0.01;
freqFactor(4) = 0.01;
freqFactor(5) = 0.01;
freqFactor(6) = 0.01;
freqFactor(7) = 0.01;
freqFactor(8) = 0.01;
freqFactor(9) = 0.01;
freqFactor(10) = 0.01;


%             
% freqFactor = [1;
%               1E-5;
%               1E-2;
%               10;
%               0.01;
%               0.01;
%              ];
%           
       
%Initialize Bat Algorithm related parameters

PopSize = 20; % no. of bats doing the search
GenSize = 200; % no. of generations of pulses.
Loudness(1:PopSize) = 0.95;
pulseRate0(1:PopSize) = 0.95;
Fmin = zeros(10,1);
Fmax = 75 * freqFactor;
d = 10; % dimensions of the solution - the parameters under scanner
% ParameterVals  = [5,5,5,3,3,3,1,1,10,10];
% ParameterVals = [6.23838337751987,4.08231143093526,2.65528420692094,3.69013962891634,3.54927049982546,3.04233374621455,1.30789351389086,1.39506660757733,14.9255287416866,9.04766383710071;];
ParameterVals = [8.41748435115720,5.11059675876186,0.384694552176183,2.33877330020370,0.354969777919614,3.01047791904419,1.32375374137403,2.39160040374453,27.6198898953542,10.3285700566350;];
% ParameterVals = [5.91688609760111,2.74631376303794,3.40720926911785,2.07742215018364,3.12157500235404,3.69633374813812,0.882811294214287,1.49012996974544,21.0463310266743,4.59473215050972;];
% Lower limit/bounds/ a vector
Lb= 0.1*ParameterVals;


% Upper limit/bounds/ a vector
Ub= 1.9*ParameterVals;

% Initializing arrays
Q=zeros(PopSize,1);   % Frequency
v=zeros(PopSize,d);   % Velocities

% Initialize the population/solutions
for i=1:PopSize,
  Sol(i,:) = Lb+(Ub-Lb).*rand(1,d);
  Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
  Fitness(i,:)= LQR_QR_Tuning(Sol(i,:));
  while (Fitness(i,1) < (1/7055)) ||  (Fitness(i,2) > 0.6) 
      Sol(i,:) = Lb+(Ub-Lb).*rand(1,d);
      Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
      Fitness(i,:)= LQR_QR_Tuning(Sol(i,:));
  end
i 
Fitness (i,:)
end

[fmin,I]=max(Fitness(:,1));

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
              S(i,:)=best+ (0.001*pulseFactor * randn(1,d)')';
          end

     % Evaluate new solutions
           S(i,:) = simplebounds(S(i,:),Lb,Ub);
           Fnew(1:2)= LQR_QR_Tuning(S(i,:));
     % Update if the solution improves, or not too loud
           if ((Fnew(1) >= Fitness(i,1))&&(Fnew(2) < 0.6)) & (rand<Loudness (i)) ,
                Sol(i,:)=S(i,:);
                Fitness(i,:)=Fnew(:);
           end

          % Update the current best solution
          if Fnew(1)>fmin && (Fnew(2) < 0.6)
                best=S(i,:);
                fmin=Fnew(1);
                FitnessStore (Fprogress) = Fnew(1);
                BestfitParam(Fprogress,:) = best(:);
                Fprogress = Fprogress + 1;
                ischanged(i) = 1;
                best
                ParameterVals
                Fnew
                
          end
        end
        
%         for i=1:PopSize,
%             Sol(i,:) = best;
%         end
        N_iter=N_iter+PopSize;
end


% save('BatTuneQRHVAC1','BestfitParam','FitnessStore');
