%% Simulate 2, generate results like Table 2. The cluster number is known.
% Clustering results for time series generated by two AR(1) models with the
%   same AR coefficient distribution range but different noise variance.

%Two AR(1) models:
% component    AR coefficient    noise variance
%     1         0.30 ± 0.01          0.01
%     2         0.30 ± 0.01          0.02

    clear;
    
    result2 = zeros(1,3); % store similarities
    G = cell(2,1); % two ground clusters
    G{1} = [1:15];
    G{2} = [16:30];
    sim = zeros(10,1);
    for t = 1:10 % 10 trils for each dataset
        phi1=0.3-0.01+(0.3+0.01-(0.3-0.01)).*rand(15,1);
        for i=1:15
            model1{i,1} = arima('Constant',0,'AR',phi1(i,1),'ARLags',1,'Variance',0.01);
            D{i,1} = simulate(model1{i,1},randi([500,1500],1));
        end
        % cluster 2
        phi2=0.3-0.01+(0.3+0.01-(0.3-0.01)).*rand(15,1);
        for i=1:15
            model2{i,1} = arima('Constant',0,'AR',phi2(i,1),'ARLags',1,'Variance',0.02);
            D{i+15,1} =simulate(model2{i,1},randi([500,1500],1));
        end
        clusterObj=cluster(D);
        clusterObj.initialize(1,0,2);
        [A,~] = clusterObj.EM();
        sim(t) = evaluate(G,A);
    end
    result2(1,1) = min(sim); % min similarity
    result2(1,2) = sum(sim)/10; % average similarity
    result2(1,3) = max(sim); % max similarity
    result2