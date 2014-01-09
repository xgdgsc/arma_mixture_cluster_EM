%% Class that implements the EM algorithm in the paper
classdef cluster < handle
    
    properties
        %class properties
        D_cell;        % time series data array, contains all the input time series data
        D_size;
        Theta_array;    % cluster parameter array
        cluster_num;
        P_w_x_Theta;  % store P_wk_xi_Theta [cluster_num * D_size]
        P_w;
    end
    
    methods
        %class methods
       %% constructor
        function obj=cluster(D) 
            obj.D_cell=D;
            [obj.D_size,~]=size(D); %get size of time series
        end
        
       %% random initialize Parameter set Theta and P arrays
        function initialize(self,p,q,known_cluster_num)
            %initialize p_w
            pw=1/known_cluster_num;
            self.P_w=ones(known_cluster_num,1).*pw;
            self.cluster_num=known_cluster_num;
            %initialize Theta_array
            self.Theta_array=cell(known_cluster_num,1);
            for i=1:known_cluster_num
                tmp_obj=parameter(p,q);
                tmp_obj.random_init();
                self.Theta_array{i,1}=tmp_obj;
            end
            %initialize P_wk_xi_Theta
            self.P_w_x_Theta=zeros([known_cluster_num,self.D_size]);
        end
        
       %% EM procedure
        function EM(self)
            
            while true
                [self.cluster_num,~]=size(self.Theta_array);
                %P_w_x_Theta=zeros(self.cluster_num,self.D_size);
                
                %% E step
                %calculate p_x_phi, prepare the calculation of p_w_x_Theta
                p_x_phi=zeros([self.D_size,self.cluster_num]);
                for i=1:self.D_size
                    for j=1:self.cluster_num
                        p_x_phi(i,j)=self.calc_P_x_phi(self.D_cell{i,1},self.Theta_array{j,1});
                    end
                end
                
                
                %Calculate P(wk|xi,Theta) as in (3) on page 4
                for i=1:self.D_size
                    denominator=0;
                    for k=1:self.cluster_num
                        self.P_w_x_Theta(k,i)=p_x_phi(i,k)*self.P_w(k,1);
                        denominator=denominator+self.P_w_x_Theta(k,i);
                    end
                    self.P_w_x_Theta(:,i)=self.P_w_x_Theta(:,i)/denominator;
                end
                
                %Calculate Q(Theta|Theta(t)) as in (2) on page 4
                Q=0;
                for i=1:self.D_size
                    for k=1:self.cluster_num
                        Q=Q+self.P_w_x_Theta(k,i)*(log(p_x_phi(i,k))+log(self.P_w(k,1)));
                    end
                end
                
                
            end
        end
        
       %% E step functions
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