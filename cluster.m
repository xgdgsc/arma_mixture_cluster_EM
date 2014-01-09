%% Class that implements the EM algorithm in the paper
classdef cluster < handle
    
    properties
        %class properties
        D_cell;        % time series data array, contains all the input time series data
        D_size;
        Theta_array;    % cluster parameter array
        P_wk_xi_Theta;  % store P_wk_xi_Theta [cluster_num * D_size(1)]
    end
    
    methods
        %class methods
        
        function obj=cluster(D) 
            %constructor
            obj.D_cell=D;
            obj.D_size=size(D); %get size of time series
        end
        
        % random initialize Parameter set Theta
        function initialize(self,p,q,known_cluster_num)
            self.Theta_array=cell(known_cluster_num,1);
            for i=1:known_cluster_num
                tmp_obj=parameter(p,q);
                tmp_obj.random_init();
                self.Theta_array{i,1}=tmp_obj;
            end
            self.P_wk_xi_Theta=zeros([known_cluster_num,self.D_size(1)]);
        end
        
        function EM(self)
            [cluster_num,~]=size(self.Theta_array);
            while true
                
%                 for k=1:cluster_num
%                     for 
%                 end
            end
        end
        
        %Calculate et (page 3, left above 3.2 )
        function et=calc_Et(self,x,t,e,parameterObj)
            sum_phi_x=0;
            for j=1:parameterObj.p
                if t-j>0
                    sum_phi_x = sum_phi_x+parameterObj.phi_array(j,1)*x(t-j,1);
                else
                    break;
                end
            end
            sum_theta_e=0;
            for j=1:parameterObj.q
                if t-j>0
                    sum_theta_e = sum_theta_e+parameterObj.theta_array(j,1)*e(t-j,1);
                else
                    break;
                end
            end
            et=x(t,1)-parameterObj.phi0-sum_phi_x-sum_theta_e;
        end
        
        %Calculate p(x|phi) (page 3 , above et=, adapted from ln())
        function p_x_phi=calc_P_x_phi(self,x,parameterObj)
           [n,~]=size(x);
           e=zeros(n,1);
           for t=1:n
               e(t,1)=self.calc_Et(x,t,e,parameterObj);
           end
           p_x_phi=exp(-n/2*log(2*pi*(parameterObj.sigma^2))-1/(2*parameterObj.sigma^2)*sum(e.^2));
        end
    end
    
end