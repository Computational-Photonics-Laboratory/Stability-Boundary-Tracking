function [theta_newpoints] = locate_check_polar(theta_points,dtheta,err,theta_dir,bif)

    %modifies the theta point set passed in based on the error and
    %bifurcation type

    if bif == 1 && (err == 1 || err == 2)
        dis = 0.25;
        if err == 1
            dis = 0.5;
        end
        theta_newpoints = theta_points - theta_dir * dis * dtheta;
        
    elseif bif == 1 && err == 3
        dis = 0.35;
        theta_newpoints = theta_points + theta_dir * dis * dtheta;
        
    elseif bif == 2
        dis = 0.35;
        if err == 2
            dis = 0.5;
        elseif err == 1
            dis = 0.75;
        end
        theta_newpoints = theta_points - theta_dir * dis * dtheta;   
    end
    
end