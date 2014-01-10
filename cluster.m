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
            % TODO:Vectorize
            while true
                [self.cluster_num,~]=size(self.Theta_array);
                %P_w_x_Theta=zeros(self.cluster_num,self.D_size);
                
              %% E step
                %calculate p_x_phi, prepare the calculation of p_w_x_Theta
                p_x_phi=zeros([self.D_size,self.cluster_num]);
                %store E cell array([D_size*cluster_num]) for M step calculation of sigma
                e_cell=cell(self.D_size,self.cluster_num);
                for i=1:self.D_size
                    for j=1:self.cluster_num
                        [e_cell{i,j},p_x_phi(i,j)]=self.calc_P_x_phi(self.D_cell{i,1},self.Theta_array{j,1});
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
                
              %% M step
                % update P_w as (A.2) on page 9
                self.P_w = mean(self.P_w_x_Theta,2);
                
                % update sigma as (A.3) on page 9
                for k=1:self.cluster_num
                    denominator=0;
                    numerator=0;
                    for i=1:self.D_size
                        [n,~]=size(self.D_cell{i,1});
                        denominator=denominator+n*self.P_w_x_Theta(k,i);
                        sum_e_2=sum(e_cell{i,k}.^2);
                        numerator=numerator+self.P_w_x_Theta(k,i)*sum_e_2;
                    end
                    sigma=sqrt(numerator/denominator);
                    self.Theta_array{k,1}.sigma=sigma;
                end
                
                % update other parameters as (A.7) on page 10
                for k=1:self.cluster_num
                    % calculate Wk
                    Wk=0;
                    for i=1:self.D_size
                        A=self.calcA(i,k);
                        B=self.calcB(e_cell,i,k);
                    end
                end
                
                
              %% Convergence condition
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
        function [e,p_x_phi]=calc_P_x_phi(self,x,parameterObj)
           [n,~]=size(x);
           e=zeros(n,1);
           for t=1:n
               e(t,1)=self.calc_Et(x,t,e,parameterObj);
           end
           p_x_phi=exp(-n/2*log(2*pi*(parameterObj.sigma^2))-1/(2*parameterObj.sigma^2)*sum(e.^2));
        end
        
        %% M step functions
        
         % calc A as shown in page 10 right column
         function A=calcA(self,i,k)
             p=self.Theta_array{k,1}.p;
             A=zeros(p+1,p+1);
             [n,~]=size(self.D_cell{i,1});
             A(1,1)=n;
             % calc a_0v
             for v=1:p
                 %s for sum
                 s=0;
                 for t=v+1:n
                     s=s+self.D_cell{i,1}(t-v,1);
                 end
                 A(1,v+1)=s;
             end
             % calc a_uv
             for v=1:p
                 for u=1:p
                     s=0;
                     for t=1:n
                         if t-u>0 && t-v>0
                             s=s+self.D_cell{i,1}(t-u,1)+self.D_cell{i,1}(t-v,1);
                         end
                     end
                     A(u+1,v+1)=s;
                 end
             end
         end
         
         % calc B as shown in page 10 right column
         function B=calcB(self,e_cell,i,k)
             p=self.Theta_array{k,1}.p;
             q=self.Theta_array{k,1}.q;
             B=zeros(p+1,q);
             e=e_cell{i,k};
             [n,~]=size(self.D_cell{i,1});
             % calc b_0v
             for v=1:q
                 s=0;
                 for t=v+1:n
                     s=s+e(t-v,1);
                 end
                 B(1,v)=s;
             end
             
             % calc b_uv
              for u=1:p
                 for v=1:q
                     s=0;
                     for t=1:n
                         if t-u>0 && t-v>0
                            s=s+e(t-u,1)*self.D_cell{i,1}(t-v,1);
                         end
                     end
                     B(u+1,v)=s;
                 end
              end
         end
         
    end
    
end