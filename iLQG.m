classdef iLQG%真的是iLQG吗？
    properties
        CostFun
        StateFun
    end
    methods
        function [xhat,uhat] = NStepiLQG(obj,x,u,N)
            xhat = x;uhat = u;
            for i = 1:N
                [xhat,uhat] = obj.OneStepiLQG(xhat,uhat);
            end
        end
        function [xhat,uhat] = OneStepiLQG(obj,x,u)
            %u和x为已知的控制和轨迹
            %返回的是修正后的轨迹
            %轨迹的长度应该是控制量长度+1
            N = length(x);
            k = containers.Map();
            %开环控制的部分k
            K = containers.Map();
            J = obj.J(x,u);
            Vx = obj.CostFun.lx(x(N,:),0);
            Vxx = obj.CostFun.lxx(x(N,:),0);
            for i = 1:N-1
                lx = obj.CostFun.lx(x(N-i,:),u(N-i,:));
                lu = obj.CostFun.lu(x(N-i,:),u(N-i,:));
                luu = obj.CostFun.luu(x(N-i,:),u(N-i,:));
                lux = obj.CostFun.lux(x(N-i,:),u(N-i,:));
                lxu = lux';
                lxx = obj.CostFun.lxx(x(N-i,:),u(N-i,:));
                Fx = obj.StateFun.Fx(x(N-i,:),u(N-i,:));
                Fu = obj.StateFun.Fu(x(N-i,:),u(N-i,:));
                Fxu = obj.StateFun.Fxu(x(N-i,:),u(N-i,:));
                Fuu = obj.StateFun.Fuu(x(N-i,:),u(N-i,:));
                Fxx = obj.StateFun.Fxx(x(N-i,:),u(N-i,:));
                Qx = lx+Fx'*Vx;
                Qu = lu+Fu'*Vx;
                Qxx = lxx+Fx'*Vxx*Fx+obj.dot(Fxx,Vx);
                Quu = luu+Fu'*Vxx*Fu+obj.dot(Fuu,Vx);
                Qxu = lxu+Fx'*Vxx*Fu+obj.dot(Fxu,Vx);
                Qux = Qxu';
                K(string(N-i)) = -inv(Quu+1e-1*diag(ones(1,length(Quu))))*Qux;
                k(string(N-i)) = -inv(Quu+1e-1*diag(ones(1,length(Quu))))*Qu;
                %更新Vxx和Vx
                Vx = Qx+Qxu*k(string(N-i));
                Vxx = Qxx+Qxu*K(string(N-i));
            end
            %更新轨迹
            xhat = x(1,:);
            %进行线搜索
            alpha = 1;
            IsBreak = false;
            for count = 1:10%差不多得了，别把电脑搞死循环了
                for i = 1:N-1
                    if(IsBreak)
                        break;
                    end
                     if(count ==10)
                         alpha = 1;
                     end
                    uhat(i,:) = u(i,:)+alpha*(k(string(i))+K(string(i))*(xhat(i,:)-x(i,:))')';
                    xhat(i+1,:) = obj.StateFun.F(xhat(i,:),uhat(i,:));
                end
                if(obj.J(xhat,uhat)>J&&count<10)
                    alpha = alpha/2;
                    %disp("当前步长是"+string(alpha))
                else 
                    IsBreak = true;
                end
            end
        end
        function A = dot(obj,F_xx,Vx)
            %返回一个矩阵，其中的每一个元素是F_xx的第三个维度
            %和Vx求和的结果
            Vx = Vx(:);
            A = sum(F_xx.*reshape(Vx, 1, 1, []),3);
        end
        function J = J(obj,x,u)
            N = length(x);
            J = obj.CostFun.l(x(N,:),0);
            for i = 1:N-1
                J = J+obj.CostFun.l(x(N-i,:),u(N-i,:)); 
            end
        end
        function obj = iLQG(CostFun_,StateFun_)
            obj.CostFun = CostFun_;
            obj.StateFun = StateFun_;
        end
    end
end
