%% parameter class to contain ARMA parameters (for one cluster)
classdef parameter < handle
    
    properties
        %class properties
        p;
        q;
        phi0;
        phi_array;
        theta_array;
        sigma;
        
    end
    
    methods
        function obj=parameter(p,q)
            obj.p=p;
            obj.q=q;
        end
        
        function random_init(self)
            self.phi0=rand;
            if self.p~=0
                self.phi_array=rand(self.p,1);
            end
            if self.q~=0
                self.theta_array=rand(self.q,1);
            end
            self.sigma=rand;
        end
    end
    
end