classdef Dynamics
       
    properties
       N;
       d;
        
  
       delta;
       
       alpha1;
     
       alpha3;
       
       alpha5;
       

    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, delta, alpha1, alpha3, alpha5)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
           
            

            obj.alpha1 = alpha1;
            obj.alpha3 = alpha3;
            obj.alpha5 = alpha5;
        end
        

        
        
        
        function res = fx(obj, v)
            res = v;
        end
        
            
        
        function res = fv(obj, x, v, u)
            res = zeros(obj.N, obj.d);

            for i=1:obj.N
                temp = zeros(1, obj.d);
                for j=1:obj.N
                    temp = temp+  obj.a(norm(x(i, :) - x(j, :))) * (v(j, :) - v(i, :));
                end 
                res(i, :) = temp/obj.N + u(i, :);
            end
        end
        
                
       
        
        function res = a(obj, r)
            res = 1/(1 + r^2)^obj.delta;
        end
        
        function res = da(obj,  r)
            res = -obj.delta*2*r / (1 + r^2)^(1 + obj.delta);
        end
        
               
        
        
        function res = fz(obj, x, v, u)
            res = obj.alpha1 * B(v, obj.N, obj.d); 
            res = res+  0.5*obj.alpha3 * norm(u)^2;
            res = res+  obj.alpha5 * B(x, obj.N, obj.d);
        end
        
        
        
        
        
        
        function res = dfvdx(obj, x, v, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+    (v(j, :) - v(i, :))' * obj.da(norm(x(i, :) - x(j, :))) * dnormdiff(x, i, j, k, obj.d);
            end
            temp = temp/obj.N;
            
            res = temp;
        end
        
       
        

        
        
        
        
        
        
        
        function res = dfvdv(obj, x, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+   obj.a(norm(x(i, :) - x(j, :))) * (dwdw(j, k, obj.d) - dwdw(i, k, obj.d));
            end
            temp = temp/obj.N;
            
            res = temp;
        end
        
                
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        function res = dfvdu(obj, i, k)
            res = zeros(obj.d, 1);
            if k == i
                res = eye(obj.d);
            end
        end
        
        function res = dfzdu(obj, u, k)
            res = obj.alpha3*u(k, :);
        end
        
        
        

    
        
        
        function res = F(obj, argx, argu)
            [x, v, z, u] = convert(argx, argu, obj.N, obj.d);
            
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(obj.fv(x, v, u)', [obj.N*obj.d, 1]);     obj.fz(x, v, u)];
        end
        
        
        function res = GuF(obj, argx, argu)
            [x, v, z, u] = convert(argx, argu, obj.N, obj.d);
            
            N = obj.N;
            d = obj.d;
            
            res = zeros(2*N*d + 1, N*d);
            
                        
%             dfvdu
            res(N*d+1:2*N*d, 1:N*d) = eye(N*d);
            
%             dfzdu
            if obj.alpha3 ~= 0
                temp = zeros(1, N*d);
                for k = 1:N
                    temp((k-1)*d+1:k*d) = obj.dfzdu(u, k); 
                end
                res(2*N*d+1, :) = temp;
            end
        end
        
        
        
        function res = GxF(obj, argx, argu)
            [x, v, z, u] = convert(argx, argu, obj.N, obj.d);
            
            N = obj.N;
            d = obj.d;
            
            temp = zeros(N*d);
            
            res = zeros(2*N*d + 1);
            
            %dfxdv
            res(1:N*d,    N*d+1:2*N*d) = eye(N*d);
                                   
            %dfvdx
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdx(x, v, i, k);
                end
            end
            res(N*d+1:2*N*d, 1:N*d) = temp;
            
            %dfvdv
            for i = 1:N
                for k=1:N
                    temp((i-1)*d+1:i*d, (k-1)*d+1:k*d) = obj.dfvdv(x, i, k);
                end
            end
            res(N*d+1:2*N*d, N*d+1:2*N*d) = temp;
            
            %dfzdv
            if obj.alpha1 ~=0
                temp = zeros(1, N*d);
                for k = 1:N
                    temp((k-1)*d+1:k*d) = obj.dfzdv(v, k); 
                end
                res(2*N*d+1, N*d+1:2*N*d) = temp;
            end
            
            %dfzdx
            if obj.alpha5 ~= 0
                temp = zeros(1, N*d);
                for k = 1:N
                    temp((k-1)*d+1:k*d) = obj.dfzdx(x, k); 
                end
                res(2*N*d+1, 1:N*d) = temp;
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        function res = dfzdv(obj, v, k)
            res = obj.alpha1 * dBdw(v, k, obj.N, obj.d);
        end
        
        function res = dfzdx(obj, x, k)
            res = obj.alpha5 * dBdw(x, k, obj.N, obj.d);
        end
        
       
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

