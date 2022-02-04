function [alpha_points,F_points,theta_points,pos1,pos2,dtheta] = find_nextboundary(alpha,F,range,r,step_num,theta_dir)

    %find boundary point in distance closest to range
    min = 1e3;
    pos1 = 1;
    pos2 = 0;
    for n = 1:length(r)
        diff = abs(range - r(n));
        if diff < min
            pos1 = n;
            min = diff;
            if range - r(n) > 0
                pos2 = n + 1;
            elseif range - r(n) < 0
                pos2 = n - 1;
            end
        end
    end

    %get initial angle guess
    alpha_diff = alpha(end-pos1)-alpha(end); F_diff = F(end-pos1)-F(end);
    theta1 = mod(atan2(F_diff,alpha_diff),2*pi);

    %repeat for second angle if applicable
    theta2 = 0;
    if pos2 ~= 0
        alpha_diff = alpha(end-pos2)-alpha(end); F_diff = F(end-pos2)-F(end);
        theta2 = mod(atan2(F_diff,alpha_diff),2*pi);
    end
    
    %number of circumference points
    num_points = step_num;
    %rad distace between points
    dtheta = 2 * pi / num_points;
    
    %array of angles of circumference points
    theta_points = zeros(2,num_points-1);
    for n = 1:num_points-1
        theta_points(1,n) = theta1 + theta_dir * dtheta * n;
        if pos2 ~= 0
            theta_points(2,n) = theta2 + theta_dir * dtheta * n;
        end
    end
    
    %convert angles to coordinates
    alpha_points = zeros(2,num_points-1);
    F_points = zeros(2,num_points-1);
    for n = 1:num_points-1
        alpha_points(1,n) = alpha(end) + range * cos(theta_points(1,n));
        F_points(1,n) = F(end) + range * sin(theta_points(1,n));
        if pos2 ~= 0
            alpha_points(2,n) = alpha(end) + range * cos(theta_points(2,n));
            F_points(2,n) = F(end) + range * sin(theta_points(2,n));
        end
    end
   
end