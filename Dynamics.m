classdef Dynamics
       
    properties
       N;
       d;
        
              
       delta;
       
       
       
       alpha1;% integral of V(t)
     
%        alpha3; % control cost
       
       alpha5;% integral of X(t)
       
       alpha7;% the potential 
       
             
       R;    
       
        
    end
    
    properties (Constant)
        Me = 1; % index of the controlled agent is always 1
    end
    
    
    
    
    
    
    methods
        
        function obj = Dynamics(N, d, delta, alpha1, alpha7, R)
            obj.N = N;
            obj.d = d;
            
            obj.delta = delta;
            
            obj.alpha1 = alpha1;% integral of Bp(v)
%             obj.alpha3 = alpha3; %control cost
            obj.alpha7 = alpha7;% the potential 
            
            obj.R = R;
            
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
                res(i, :) = temp/obj.N; 
            end
            
            res(obj.Me, :) = res(obj.Me, :)+    u;
        end
        
                
       
        
        function res = a(obj, r)
            res = 1/(1 + r^2)^obj.delta;
        end
        
        function res = da(obj,  r)
            res = -obj.delta*2*r / (1 + r^2)^(1 + obj.delta);
        end
        
               
        
        
        function res = fz(obj, x, v)
            res = obj.alpha1 * Bp(v, obj.Me, obj.N); 
            
            
            if obj.alpha7 ~= 0
                temp = 0;
                for j = 1:obj.N 
                    if j ~= obj.Me
                        temp = temp+    V(norm(x(obj.Me, :) - x(j, :)), obj.R, obj.N);
                    end
                end
                res = res+  obj.alpha7 * temp;
            end
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
%             dfvdu            
            res = eye(obj.d, obj.d);
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
            if obj.alpha7 ~= 0
                temp = zeros(1, N*d);
                for k = 1:N
                    temp((k-1)*d+1:k*d) = obj.dfzdx(x, k); 
                end
                res(2*N*d+1, 1:N*d) = temp;
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        function res = dfzdv(obj, v, k)
            res = obj.alpha1 * dBpdw(v, obj.Me, k, obj.N, obj.d);
        end
        
        function res = dfzdx(obj, x, k)
            res = zeros(1, obj.d);
            
            if obj.alpha7 ~= 0
                temp = zeros(1, obj.d);
                for j = 1:obj.N
                    if j ~= obj.Me
                        temp = temp+    dV(norm(x(obj.Me, :) - x(j, :)), obj.R, obj.N) * dnormdiff(x, obj.Me, j, k, obj.d);
                    end
                end
                res = res+  obj.alpha7 *temp;
            end
        end
               
    
        
        
    end 
    
    
    methods(Static)
        
        
    end
    
end

