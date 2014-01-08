%% Class that implements the EM algorithm in the paper
classdef cluster 
    properties
        D_array;% time series data array, contains all the input time series data
    end
    
    methods
        function obj=cluster(D)
            obj.D_array=D;
        end
    end
    
end