function [u] = convert_control(argu, Nu, d)
            u = reshape(argu, [d, Nu])';
end
