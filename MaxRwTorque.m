function [r, Tmax_L2] = MaxRwTorque(v, W, C)
% This function finds the maximum torque in 'v' direction of the RWA whose
% torque distribution matrix is given by W.
%
% Written by Hyosang Yoon (hyosang.yoon@kaist.ac.kr)
% at KAIST, Dept. of Aerospace Engineering
% 
% Copyrightⓒ2020. Hyosang Yoon. All rights reserved.
%
% INPUT
% v: torque direction
% W: RW torque distribution matrix
% C: coplanar RW combinations, output of 'CplnrWhlsCombn' function
%
% OUTPUT
% r: RW torque ratio, should be in [-1, 1] for causal systems
% Tmax_L2: magnitude (L2 norm) of the maximum RW torque in 'v' direction

machine_zero = 1e-10; 
v = v / norm(v);
for i = 1:1:length(C)
    a = C{i}(1);
    b = C{i}(2);
    
    ni = cross(W(:,a),W(:,b));              % Eq.(6)
    ni = ni/norm(ni);           
    if(dot(ni, v) > machine_zero)           % Eq.(8)
        r = sign(W'*ni);
    elseif(dot(ni, v) < -machine_zero)
        r = -sign(W'*ni);
    else
        % If the plane is parallel to the torque direction, it cannot be 
        % the bounding facet.
        continue;
    end
    r(C{i}) = zeros(length(C{i}),1); 
    
    Tk = W * r;                             % Eq.(9)
    Tmax_L2 = dot(ni, Tk) / dot(ni, v);     % Eq.(10)-(12)
    Tmax = Tmax_L2 * v;                     % Eq.(13)
    Tc = Tmax - Tk;                         % Eq.(14)
    if(norm(Tc) < machine_zero)
        % If Tc = 0, rc should be zero
        return;
    end
    
    Wc = W(:,C{i});
    [rc] = find_rc(Tc, Wc, ni);
    if(all(abs(rc) < 1 + machine_zero))
        r(C{i}) = rc;
        return
    end
end

end


function [rc] = find_rc(Tc, Wc, ni)
machine_zero = 1e-10;
N = size(Wc, 2);
Tc_L2 = norm(Tc);
vc = Tc / Tc_L2;                                % Eq.(17)
for j = 1:1:N
    nj = cross(Wc(:,j), ni) / norm(Wc(:,j));    % Eq.(18)
    if(dot(nj, vc) > machine_zero)              % Eq.(16)
        rcmax = sign(Wc' * nj); 
    elseif(dot(nj, vc) < -machine_zero)
        rcmax = -sign(Wc' * nj);
    else
        % If the edge is parallel to the torque direction, it cannot be 
        % the bounding edge.
        continue;
    end
    rcmax(j) = 0;

    TL = Wc * rcmax;                                    % Eq.(19)
    Tcmax_L2 = dot(nj, TL) / dot(nj, vc);               % Eq.(20)-(22)
    Tcmax = Tcmax_L2 * vc;                              % Eq.(23)
    TJ = Tcmax - TL;                                    % Eq.(24)
    rcmax(j) = dot(TJ, Wc(:,j)) / dot(Wc(:,j),Wc(:,j)); % Eq.(25)
    rc = Tc_L2 / Tcmax_L2 * rcmax;                      % Eq.(26)
    if(all(abs(rc) < 1 + machine_zero))
        return
    end
end

end
