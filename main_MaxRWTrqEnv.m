%% 
% Example code for "Maximum Reaction-wheel Array Torque/Momentum Envelopes 
% for General Configurations"
%
% Written by Hyosang Yoon (hyosang.yoon@kaist.ac.kr)
% at KAIST, Dept. of Aerospace Engineering
%
% Copyrightⓒ2020. Hyosang Yoon. All rights reserved.

%% Reaction wheel configurations
clear
clc
close all

cases = cell(1,3);

cases{1}.A = [...
    1, -1, 1, -1;
    1, 1, -1, -1;
    1, 1, 1, 1];
cases{1}.Tmax = [100,100,100,100];

cases{2}.A = [...
    1, -1, 1, -1, 1;
    1, 1, -1, -1, 0;
    1, 1, 1, 1, 0];
cases{2}.Tmax = [100,100,100,100,100];

cases{3}.A = [...
    1, -1, 1, -1, 1, 0, 1, 1;
    1, 1, -1, -1, 0, 1, 1, -1;
    1, 1, 1, 1, 0, 0, 0, 0];
cases{3}.Tmax = [10,20,30,40,50,60,70,80];

for i = 1:1:length(cases)
    A = cases{i}.A;
    RWN = size(A, 2);
    for j = 1:1:RWN
        A(:,j) = A(:,j) / norm(A(:,j));
    end
    W = A * diag(cases{i}.Tmax);
    cases{i}.W = W;
end

%% Draw the torque envelopes
for c = 1:1:length(cases)
    W = cases{c}.W;
    RWN = size(W, 2);
    [C, K] = CplnrWhlsCombn(W);

    [X,Y,Z] = sphere(64);
    N = size(X,1);

    for i = 1:1:N
        for j = 1:1:N
            v = [X(i,j);Y(i,j);Z(i,j)];
            [r, Tmax_L2] = MaxRwTorque(v, W, C);
            v = W * r;
            X(i,j) = v(1);
            Y(i,j) = v(2);
            Z(i,j) = v(3);
        end
    end
    
    figure
    surf(X,Y,Z, 'EdgeColor', 'k', 'FaceColor', 'w')
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title(['Torque Envelope for ', num2str(RWN), ' RWs'])

end