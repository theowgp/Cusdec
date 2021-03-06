classdef Objective
    
    
    properties
        N;
        d;
        
        alpha2;
       
        alpha4;
        
        dynamics;
    end
    
    methods
        function obj = Objective(dynamics, N, d,  alpha2,  alpha4)
            obj.N = N;
            obj.d = d;
            
            
            obj.alpha2 = alpha2;
            obj.alpha4 = alpha4;
            
            obj.dynamics = dynamics;
        end
        
        
%         function res = phi(obj, z)
%             res = z;
%         end
        function res = phi(obj, argx)
            [x, v, z] = convert_state(argx, obj.N, obj.d);
            res = z;
            res = res+  B(v, v, obj.N) * obj.alpha2;
            res = res+  E(x, v, obj.N) * obj.alpha4;
        end
         
        function res = Gxphi(obj, argx)
            [x, v, z] = convert_state(argx, obj.N, obj.d);
            res = zeros(2*obj.N*obj.d + 1, 1);
            
%             dphidx
            if obj.alpha4 ~= 0
                temp1 = zeros(obj.N*obj.d, 1);
                for k = 1:obj.N
                    temp1((k-1)*obj.d+1:k*obj.d, 1) = dBdw(x, k, obj.N, obj.d); 
                end

                X = B(x, x, obj.N);

                res(1:obj.N*obj.d, 1) = obj.dynamics.a(sqrt(2*obj.N*X)) * 0.5 * temp1 / sqrt(X) * obj.alpha4;
            end
            
%             ddhidv
            if obj.alpha2 ~= 0  || obj.alpha4 ~= 0
                temp2 = zeros(obj.N*obj.d, 1);
                for k = 1:obj.N
                    temp2((k-1)*obj.d+1:k*obj.d, 1) = dBdw(v, k, obj.N, obj.d); 
                end

                V = B(v, v, obj.N);

                res(obj.N*obj.d+1:2*obj.N*obj.d, 1) = temp2*obj.alpha2    +   0.5 * temp2 / sqrt(V) * obj.alpha4; 
            end
            
            
%             dphidz
            res(2*obj.N*obj.d + 1) = 1;
        end
    end
    
end
