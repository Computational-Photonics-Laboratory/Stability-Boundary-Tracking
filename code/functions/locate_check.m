function [F_newpoints] = locate_check(dF,Fdir,err,F_points,bif)
    
    %modifies the F point set passed in based on the error and bifurcation
    %type
    
    if bif == 1 && (err == 1 || err == 2)
        dis = 0.25;
        if err == 1
            dis = 0.5;
        end
        F_newpoints = F_points - Fdir * dis * dF;
        
    elseif bif == 1 && err == 3
        dis = 0.35;
        F_newpoints = F_points + Fdir * dis * dF;
        
    elseif bif == 2
        dis = 0.25;
        if err == 2
            dis = 0.5;
        elseif err == 1
            dis = 0.75;
        end
        F_newpoints = F_points - Fdir * dis * dF;
    end
    
end