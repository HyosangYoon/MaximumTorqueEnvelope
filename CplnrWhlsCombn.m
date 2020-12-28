function [C, K] = CplnrWhlsCombn(W)
% This function finds coplanar-wheel combinations
%
% Written by Hyosang Yoon (hyosang.yoon@kaist.ac.kr)
% at KAIST, Dept. of Aerospace Engineering
%
% Copyrightⓒ2020. Hyosang Yoon. All rights reserved.
%
% Input
% W: RW torque distribution matrix
%
% Output
% C: coplanar RW combination (C_i)
% K: U - C_i

RWN = size(W, 2);
machine_zero = 1e-10;
check_matrix = zeros(RWN*(RWN-1)/2, RWN);
C = cell(1,RWN * (RWN-1)/2);
K = cell(1,RWN * (RWN-1)/2);
ii = 0;
for check_rwn = (RWN-1):-1:2
    combos = nchoosek(1:RWN,check_rwn);
    Ncombos = size(combos, 1);
    for i = 1:1:Ncombos
        % Check if the combination is included in the previous Ci-s already
        passthiscomb = 0;
        for j = 1:1:ii
            same_i = 0;
            for k = 1:1:check_rwn
                if(check_matrix(j, combos(i,k)))
                    same_i = same_i + 1;
                else
                    break;
                end
            end
            if(same_i == check_rwn)
                passthiscomb = 1;
                break;
            end
        end
        if(passthiscomb)
            continue;
        end
        
        % Check if the RWs in the combination are coplanar
        nvec = cross(W(:,combos(i,1)), W(:,combos(i,2)));
        nvec = nvec / norm(nvec);
        is_cp = 1;
        for j = 3:1:check_rwn
            if(abs(dot(nvec, W(:,combos(i,j)))) > machine_zero)
                is_cp = 0;
                break;
            end
        end
        if(is_cp)
            ii = ii + 1;
            check_matrix(ii, combos(i,:)) = 1;
            a = 1:RWN;
            a(combos(i,:)) = 0;
            
            % Ci-s and Ki-s
            C{ii} = combos(i,:);
            K{ii} = find(a);
        end
    end
end

C = C(1:ii);
K = K(1:ii);

end