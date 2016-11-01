classdef Objective
    
    
    properties
        N;
        d;
        
        alpha2;% V(T)
       
        
        
        dynamics;
    end
    
    properties (Constant)
        Me = 1; % index of the controlled agent is always 1
    end
    
    methods
        function obj = Objective(dynamics, N, d,  alpha2)
            obj.N = N;
            obj.d = d;
            
            
            obj.alpha2 = alpha2;% V(T)
            
            
            obj.dynamics = dynamics;
        end
        
        

        function res = phi(obj, argx)
            [x, v, z] = convert_state(argx, obj.N, obj.d);
            res = z;
            res = res+  Bp(v, obj.Me, obj.N) * obj.alpha2;
        end
         
        function res = Gxphi(obj, argx)
            [x, v, z] = convert_state(argx, obj.N, obj.d);
            res = zeros(2*obj.N*obj.d + 1, 1);
            
            
%             ddhidv
            if obj.alpha2 ~= 0 
                temp = zeros(obj.N*obj.d, 1);
                for k = 1:obj.N
                    temp((k-1)*obj.d+1:k*obj.d, 1) = dBpdw(v, obj.Me, k, obj.N, obj.d);
                end

                res(obj.N*obj.d+1:2*obj.N*obj.d, 1) = temp*obj.alpha2;
            end
            
            
%             dphidz
            res(2*obj.N*obj.d + 1) = 1;
        end
    end
    
end
