function [u] = convert_control(argu, N, d)
            u = reshape(argu, [d, N])';
end
