dt = 0.01;
m = 1; % 定义m的值
l = 1; % 定义l的值
g = 9.81; % 重力加速度
M = 10; % 定义M的值
J = m*l^2/3; % 定义J的值
syms x0 x1 x2 x3 u
F_ = [x0 + dt*x1;
    x1 + dt*((J + m*l^2) * u + m*l*(J + m*l^2) * sin(x2) * x3^2 - m*l^2*g ...
    *sin(x2)*cos(x2)) / ((J + m*l^2) * (M + m) - m^2*l^2*cos(x2)^2);
    x2 + dt*x3;
    x3 + dt*((M + m) * m*l*g*sin(x2) - m*l*cos(x2)*u - m^2*l^2*sin(x2)*cos(x2) * ...
    x3^2) / ((J + m*l^2) * (M + m) - m^2*l^2*cos(x2)^2);];
l_ = 10*x0^2+1e-256*x1^2+1000*x2^2+0.001*u^2+1e-15*x3^2;
StateFun = StateFunction(F_,dt);
CostFun = CostFunction(l_);
ILQG = iLQG(CostFun,StateFun);
%获取原始轨迹
N = 30;Initial_X = [0,0,pi/10,0];
Initial_U = zeros(N,1);
for i = 1:N
    Initial_X(i+1,:) = StateFun.F(Initial_X(i,:),Initial_U(i,:));
end
plot(Initial_X(:,3))
[xhat,uhat] = ILQG.NStepiLQG(Initial_X,Initial_U,20);
hold on
plot(xhat(:,3))


