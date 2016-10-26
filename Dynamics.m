classdef Dynamics
       
    properties
       N;
       d;
        
              
       delta;
       
       eps;
       
       R;
        
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, delta)
            obj.N = N;
            obj.d = d;
            obj.delta = delta;
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
        
               
        
        
        
        function res = fz(obj, v)
            res = B(v, v, obj.N);
        end
        
        
        
        
        
        
        function res = dfvdx(obj, x, v, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+    (v(j, :) - v(i, :))' * obj.da(norm(x(i, :) - x(j, :))) * obj.dnorm(x, i, j, k);
            end
            temp = temp/obj.N;
            
            res = temp;  
        end
        
        
        function res = dnorm(obj, x, i, j, k)
            res = zeros(1, obj.d);
            if i ~= j
                if k == i
                    res = x(k, :)/norm(x(k, :) - x(j, :));
                else
                    if k == j
                        res = -x(k, :)/norm(x(i, :) - x(k, :));
                    end
                end
            end
        end
        
        
        
        
        
        function res = dfvdv(obj, x, i, k)
            temp = zeros(obj.d, obj.d);
            
            for j = 1:obj.N
                temp = temp+   obj.a(norm(x(i, :) - x(j, :))) * (obj.dwdw(j, k) - obj.dwdw(i, k));
            end
            temp = temp/obj.N;
            
            res = temp;
        end
        
        
        function res = dwdw(obj, i, k)
            res = zeros(obj.d);
            if(k == i)
                res = eye(obj.d);
            end
        end
        
        
        
        
        
        
        

    
        
        
        function res = F(obj, argx, argu)
            [x, v, z, u] = convert(argx, argu, obj.N, obj.d);
            
            res = [reshape(obj.fx(v)', [obj.N*obj.d, 1]);    reshape(obj.fv(x, v, u)', [obj.N*obj.d, 1]);     obj.fz(v)];
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
            temp = zeros(1, N*d);
            for k = 1:N
                temp((k-1)*d+1:k*d) = obj.dfzdv(v, k); 
            end
            res(2*N*d+1, N*d+1:2*N*d) = temp;
        end
        
        
        
        
        function res = GuF(obj, argx, argu)
            N = obj.N;
            d = obj.d;
            
            res = zeros(2*N*d + 1, N*d);
                        
            %dfvdu
            res(N*d+1:2*N*d, 1:N*d) = eye(N*d);
        end
        
        
        
        
        
        
        
        function res = dfzdv(obj, v, k)
            res = dBdw(v, k, obj.N, obj.d);
        end
        
        
        
        
        
        
          
        
        
        
                
        
        
        
        
        
        function res = a(obj, r)
            res = 1/(1 + r^2)^obj.delta;
        end
        
        function res = da(obj,  r)
            res = -obj.delta*2*r / (1 + r^2)^(1 + obj.delta);
        end
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

