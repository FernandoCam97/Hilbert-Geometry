clear all; close all;
%% Setup
n = 100;
final_time = 7;
p_0 = [0.7;0.2;0.1];
lambda_1 = 1/2;
lambda_2 = 1;
lambda_3 = 2;
A = zeros(3,3);
A(1,1) = lambda_1;
A(2,2) = lambda_2;
A(3,3) = lambda_3;
%% Make Triangle
t = linspace(0,final_time,n);
L1 = getLine([1;0;0],[0;1;0]);
L2 = getLine([1;0;0],[0;0;1]);
L3 = getLine([0;1;0],[0;0;1]);
figure(1)
plot3(L1(1,:),L1(2,:),L1(3,:),'k','LineWidth',2)
hold on
plot3(L2(1,:),L2(2,:),L2(3,:),'k','LineWidth',2)
hold on
plot3(L3(1,:),L3(2,:),L3(3,:),'k','LineWidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
grid on
%% Build Trajectory
k = 100; % Points in tangents
Trajectory = zeros(3,n);
Trajectory2 = zeros(3,n);
Trajectory_der = zeros(3,n);
Tangents = zeros(3,k,n);
xi_m = zeros(3,n); % Intersections
xi_p = zeros(3,n);
for i = 1 : n
    Trajectory(:,i) = getTrajectory(t(i),A,p_0); 
    Trajectory2(:,i) = getTrajectory2(t(i),A,p_0); 
    % Project the point back to the triangle
    Trajectory(:,i) = Trajectory(:,i)/(getProjection(t(i),p_0,A));
    Trajectory2(:,i) = Trajectory2(:,i)/(getProjection2(t(i),p_0,A));
    Trajectory_der(:,i) = getDerivative(t(i),A,p_0);
    Tangents(:,:,i) = getLineF(Trajectory(:,i),Trajectory_der(:,i),k);
    %xi_m(:,i) = findIntersection(L1,Tangents(:,:,i),n);
    %xi_p(:,i) = findIntersection(L3,Tangents(:,:,i),n);
    
     figure(1)
     hold on
     plot3(Trajectory(1,1:i),Trajectory(2,1:i),Trajectory(3,1:i),'m')
     hold on
     plot3(Tangents(1,:,i),Tangents(2,:,i),Tangents(3,:,i),'c')
     hold on
     plot3(Trajectory2(1,1:i),Trajectory2(2,1:i),Trajectory2(3,1:i),'m')
     %scatter3(xi_m(1,1:i),xi_m(2,1:i),xi_m(3,1:i),'b')
     hold on
     %scatter3(xi_p(1,1:i),xi_p(2,1:i),xi_p(3,1:i),'b')
     pause(0.1)
end

%% Functions
function p_t = getTrajectory(t,A,p_0)
e_At = expm(A*t);
p_t = e_At*p_0;
end

function p_t = getTrajectory2(t,A,p_0)
p_t = zeros(3,1);
for i = 1 : length(A)
    p_t(i) = (A(i,i)^t)*p_0(i);
end
end

function L = getLine(p1,p2)
L = zeros(3,100);
t = linspace(0,1,100);
for i = 1 : length(t)
    L(:,i) = (1-t(i))*p1 + t(i)*p2;
end
end

function s = getProjection(t,p_0,A)
s = 0;
for i = 1 : length(p_0)
    s = s + exp(A(i,i)*t)*p_0(i);
end
end

function s = getProjection2(t,p_0,A)
s = 0;
for i = 1 : length(p_0)
    s = s + (A(i,i)^t)*p_0(i);
end
end

function s = getProjectionDer(t,p_0,A)
s = 0;
for i = 1 : length(p_0)
    s = s + A(i,i)*exp(A(i,i)*t)*p_0(i);
end
end

function p_dot = getDerivative(t,A,p_0)
p_dot = zeros(3,1);
p_dot = -(getProjection(t,p_0,A)^(-2))*(getProjectionDer(t,p_0,A))...
    *expm(A*t)*p_0 + (getProjection(t,p_0,A)^(-1))*...
    [A(1,1)*exp(A(1,1)*t)*p_0(1);A(2,2)*exp(A(2,2)*t)*p_0(2);...
    A(3,3)*exp(A(3,3)*t)*p_0(3)];
end

function L = getLineF(p,d,k)
L = zeros(3,k);
t = linspace(-1,1,k);
d = d/getNorm(d);
for i = 1 : length(L)
    L(:,i) = p + t(i)*d;
end
end

function n = getNorm(x)
n = (x(1)^2+x(2)^2+x(3)^2)^(1/2);
end

function x = findIntersection(L1,L2,n)
tol = 2/n;
for i = 1 : length(L1)
    for j = 1 : length(L2)
        if L1(1,i)-tol <= L2(1,j) && L1(1,i) + tol >= L2(1,j)
            if L1(2,i)-tol <= L2(2,j) && L1(2,i) + tol >= L2(2,j)
                if L1(3,i)-tol <= L2(3,j) && L1(3,i) + tol >= L2(3,j)
                    x = L1(:,i);
                end
            end
        end
    end
end
end

