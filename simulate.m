%% AR generating example, 
% TODO: write function wrapper to generate all the required time series
% model1 = arima('Constant',0,'AR',0.3,'ARLags',1,'Variance',0.01);
% [Y1,E]=simulate(model1,1000);
% figure(1);
% plot(Y1);
% 
% model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
% [Y2,E]=simulate(model2,1000);
% figure(2);
% plot(Y2);
clear;
%% Simulate Process
%generate test time series 

% cluster1
phi1=0.3-0.01+(0.3+0.01-(0.3-0.01)).*rand(15,1);
for i=1:15
    model1{i,1} = arima('Constant',0,'AR',phi1(i,1),'ARLags',1,'Variance',0.01);
    D{i,1} = simulate(model1{i,1},randi([500,1500],1));
end
% cluster 2
phi2=0.6-0.01+(0.6+0.01-(0.6-0.01)).*rand(15,1);
for i=1:15
    model2{i,1} = arima('Constant',0,'AR',phi2(i,1),'ARLags',1,'Variance',0.01);
    D{i+15,1} =simulate(model2{i,1},randi([500,1500],1));
end
clusterObj=cluster(D);
clusterObj.initialize(1,0,2);
clusterObj.EM();