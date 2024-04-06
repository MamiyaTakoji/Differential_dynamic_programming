classdef CostFunction
    properties
        l;
        lx;
        lu;
        luu;
        lux;
        lxx;
        %J=\frac{1}{2}x_{N}^{T}Sx_N
        % +\frac{1}{2}\sum_{k=1}^{N-1}{\left( x_{k}^{T}Qx_k+u_{k}^{T}Ru_k \right)}
    end
    methods
        function obj = CostFunction(l_)%l_为符号函数，变量为x1,x2,···,u1,u2,···
            % 使用正则表达式分离变量
            Vars = symvar(l_);
            vars = regexp(char(symvar(l_)), '\w+', 'match');
            % 分离出 x 和 u 的变量，并将其转换为 sym 类型
            x_vars = sym(vars(startsWith(vars, 'x')));
            u_vars = sym(vars(startsWith(vars, 'u')));
            obj.l = matlabFunction(l_,"Vars",{x_vars,u_vars});
            lx_ = sym([]);
            for i = 1: length(x_vars)
                lx_(end+1) = diff(l_,x_vars(i));
            end
            obj.lx = matlabFunction(lx_(:),"Vars",{x_vars,u_vars});
            lu_ = sym([]);
            for i = 1: length(u_vars)
                lu_(end+1) = diff(l_,u_vars(i));
            end
            obj.lu = matlabFunction(lu_(:),"Vars",{x_vars,u_vars});
            lxx_ = sym([]);
            for i = 1:length(x_vars)
                for j = 1:length(x_vars)
                    lxx_(i,j) = diff(diff(l_,x_vars(i)),x_vars(j));
                end
            end
            obj.lxx = matlabFunction(lxx_,"Vars",{x_vars,u_vars});
            lux_ = sym([]);
            for i = 1:length(u_vars)
                for j = 1:length(x_vars)
                    lux_(i,j) = diff(diff(l_, u_vars(i)), x_vars(j));
                end
            end
            obj.lux = matlabFunction(lux_,"Vars",{x_vars,u_vars});
            luu_ = sym([]);
            for i = 1:length(u_vars)
                for j = 1:length(u_vars)
                    luu_(i,j) = diff(diff(l_, u_vars(i)), u_vars(j));
                end
            end
            obj.luu = matlabFunction(luu_,"Vars",{x_vars,u_vars});
        end
    end
end