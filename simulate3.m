%% Simulate 3, the cluster number is unknown.
%   Generate BIC values for a typica 3-cluster time series dataset.
%   The results are like Table 4 in the paper.

%Typical 3 AR(1) models:
% component    AR coefficient    noise variance
%     1         0.20 ± 0.01      0.01 ± 0.001
%     2         0.50 ± 0.01      0.01 ± 0.001
%     3         0.80 ± 0.01      0.01 ± 0.001

clear;
maxCluster_num = 7;  % the possible maximum cluster number
sim_result = []; % store similarity under a certain cluster_num
bic_result = {}; % store bic values under a certain cluster_num

G = cell(3,1); % three ground clusters
G{1} = [1:12];
G{2} = [13:24];
G{3} = [25:36];

% model 1
phi1 = 0.19 + 0.02.*rand(12,1);
var1 = 0.009 + 0.002.*rand(12,1);
for i=1:12
    model1{i,1} = arima('Constant',0,'AR',phi1(i,1),'ARLags',1,'Variance',var1(i,1));
    D{i,1} = simulate(model1{i,1},randi([500,1500],1));
end

% model 2
phi2 = 0.49 + 0.02.*rand(12,1);
var2 = 0.009 + 0.002.*rand(12,1);
for i=1:12
    model2{i,1} = arima('Constant',0,'AR',phi1(i,1),'ARLags',1,'Variance',var2(i,1));
    D{i+12,1} = simulate(model2{i,1},randi([500,1500],1));
end

% model 3
phi3 = 0.79 + 0.02.*rand(12,1);
var3 = 0.009 + 0.002.*rand(12,1);
for i=1:12
    model3{i,1} = arima('Constant',0,'AR',phi1(i,1),'ARLags',1,'Variance',var3(i,1));
    D{i+24,1} = simulate(model3{i,1},randi([500,1500],1));
end
  
for n = 2:maxCluster_num
    clusterObj=cluster(D);
    clusterObj.initialize(1,0,n);
    [A,bic] = clusterObj.EM();
    sim_result(end+1) = evaluate(G,A);
    bic_result{end+1} = bic;
end
