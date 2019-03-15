%This file is a function where for a given set of parameters, a FIT % is
%calculated and sent bact to the PSO tuning algorithm

function fit = GetFitNL(xx)

global  B1North_R_T B1N_RWT B1N_SWT B1South_R_T B1S_RWT B1S_SWT B1Amb_TMAX    ParameterVals   
global  theta_swn theta_sws theta_o  State1 State2 State3  State4
%model Parameters
theta_swn = B1N_SWT(1,1);
theta_sws = B1S_SWT(1,1);
theta_o = B1Amb_TMAX (1,1);
ExptLength = 300;

x_ss(1:4) = [B1N_RWT(1,1); B1S_RWT(1,1);B1North_R_T(1,1);B1South_R_T(1,1)];
% x_ss(5:10) = xx(1:6);
% [t,x] = ode45('Tuning_ODE',[0 1000],x_ss);
% [i,j] = size(t);
% x_ss(1:4) = x(i,1:4);
% 
% if ((( x_ss(1) > B1N_RWT(1,1) - 5)) && ( x_ss(1) < B1N_RWT(1,1) + 5) ) ...
%     && ((( x_ss(2) > B1S_RWT(1,1) - 5)) && ( x_ss(2) < B1S_RWT(1,1) + 5) )...
%     &&((( x_ss(3) > B1North_R_T(1,1) - 5)) && ( x_ss(3) < B1North_R_T(1,1) + 5) )...
%     &&((( x_ss(4) > B1South_R_T(1,1) - 5)) && ( x_ss(4) < B1South_R_T(1,1) + 5) )
%     goAhead = 1;
% else
%     goAhead = 0;
% end


% if (goAhead == 1)
x_ss(5:10) = xx(1:6);

State1 = zeros(ExptLength,1);
State2 = zeros(ExptLength,1);
State3 = zeros(ExptLength,1);
State4 = zeros(ExptLength,1);
count = 1;

while (count <= ExptLength)
theta_swn = B1N_SWT(1,count);
theta_sws = B1S_SWT(1,count);
theta_o = B1Amb_TMAX (1,count);
 
[t,x] = ode45('Tuning_ODE',[0 18],x_ss);
[i,j] = size(t);
x_ss(1:4) = x(i,1:4);
x_ss(5:10) = xx(1:6);

State1(count) = x_ss(1);
State2(count) =  x_ss(2);
State3(count) =  x_ss(3);
State4(count) =  x_ss(4);


 count = count + 1;
end

SpredErr1 = 0;
SpredErr2 = 0;
SpredErr3 = 0;
SpredErr4 = 0;
Svariance1= 0;
Svariance2= 0;
Svariance3= 0;
Svariance4= 0;
meanS1 = mean(B1N_RWT(1:ExptLength)',1);
meanS2 = mean(B1S_RWT(1:ExptLength)',1);
meanS3 = mean(B1North_R_T(1:ExptLength)',1);
meanS4 = mean(B1South_R_T(1:ExptLength)',1);

for count = 1:ExptLength
    SpredErr1 = SpredErr1 + abs(B1N_RWT(count) - State1(count)); 
    Svariance1 = Svariance1 + abs(B1N_RWT(count) - meanS1);
    
    SpredErr2 = SpredErr2 + abs(B1S_RWT(count) - State2(count)); 
    Svariance2 = Svariance2 + abs(B1S_RWT(count) - meanS2) ;
    
    SpredErr3 = SpredErr3 + abs(B1North_R_T(count) - State3(count)); 
    Svariance3 = Svariance3 + abs(B1North_R_T(count) - meanS3); 
    
    SpredErr4 = SpredErr4 + abs(B1South_R_T(count) - State4(count)); 
    Svariance4 = Svariance4 + abs(B1South_R_T(count) - meanS4); 
    
   
end

fitVal = 100 * (1-(0.2*(SpredErr1/Svariance1)+ 0.2*(SpredErr2/Svariance2)+ 0.3*(SpredErr3/Svariance3)+ 0.3*(SpredErr4/Svariance4)));
% else
%  fitVal = 0;   
% end

% plot(Gnew,'k');
% hold on;
% plot(Gnewresp,'b');
% plot(GluReleaseresp);
fit = fitVal; 
