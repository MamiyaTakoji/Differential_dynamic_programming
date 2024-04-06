classdef StateFunction
    properties
        dt
        F
        Fx
        Fu
        Fxx
        Fxu
        Fuu
    end
    methods
        function x=OnestepIteration(obj,x_,u_)
            x(1) = x_(1)+obj.dt*x_(2);
            x(2) = x_(2)+obj.dt*(u_+sin(x_(1)));
        end
        %第一个分量是f的维数，这里为2，下同
        function obj = StateFunction(F_,dt_)
            obj.dt = dt_;
            % 使用正则表达式分离变量
            Vars = symvar(F_);
            vars = regexp(char(symvar(F_)), '\w+', 'match');
            % 分离出 x 和 u 的变量，并将其转换为 sym 类型
            x_vars = sym(vars(startsWith(vars, 'x')));
            u_vars = sym(vars(startsWith(vars, 'u')));
            obj.F = matlabFunction(F_,"Vars",{x_vars,u_vars});
            Fx_ = sym([]);
            for i = 1:length(x_vars)
                for j = 1:length(F_)
                    Fx_(j,i) = diff(F_(j),x_vars(i));
                end
            end
            obj.Fx = matlabFunction(Fx_,"Vars",{x_vars,u_vars});
            Fxu_ = sym([]);
            for i = 1:1:length(F_)
                for j = 1:length(x_vars)
                    for k = 1:length(u_vars)
                        Fxu_(j,k,i) = diff(diff(F_(i),x_vars(j)),u_vars(k));
                    end
                end
            end
            obj.Fxu = matlabFunction(Fxu_,"Vars",{x_vars,u_vars});
            Fxx_ = sym([]);
            for i = 1:1:length(F_)
                for j = 1:length(x_vars)
                    for k = 1:length(x_vars)
                        Temp = diff(F_(i),x_vars(j));
                        Fxx_(j,k,i) = diff(Temp,x_vars(k));
                    end
                end
            end
            obj.Fxx = matlabFunction(Fxx_,"Vars",{x_vars,u_vars});
            Fuu_ = sym([]);
            for i = 1:1:length(F_)
                for j = 1:length(u_vars)
                    for k = 1:length(u_vars)
                        Fuu_(j,k,i) = diff(diff(F_(i),u_vars(j)),u_vars(k));
                    end
                end
            end
           obj.Fuu = matlabFunction(Fuu_,"Vars",{x_vars,u_vars});
            Fu_ = sym([]);
            for i = 1:length(u_vars)
                for j = 1:length(F_)
                    Fu_(j,i) = diff(F_(j),u_vars(i));
                end
            end
            %[Fu_H,Fu_W] = size(Fu_);
            obj.Fu = matlabFunction(Fu_,"Vars",{x_vars,u_vars});
        end
    end
end