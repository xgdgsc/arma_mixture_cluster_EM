%% Class that implements the EM algorithm in the paper
classdef cluster < handle
    
    properties
        %class properties
        D_cell;        % time series data array, contains all the input time series data
        D_size;
        Theta_array;    % cluster parameter array
        
    end
    
    methods
        %class methods
        
        function obj=cluster(D) 
            %constructor
            obj.D_cell=D;
            obj.D_size=size(D); %get size of time series
        end
        
        %random initialize Parameter set Theta
        function initialize(self,p,q,known_cluster_num)
            self.Theta_array=cell(known_cluster_num,1);
            for i=1:known_cluster_num
                tmp_obj=parameter(p,q);
                tmp_obj.random_init();
                self.Theta_array{i,1}=tmp_obj;
            end
        end
        
        function EM(self)
        end
        
    end
    
end