function [e_val,e_sol,err] = check_eigenvalue(alpha,F_points,L,N,sol,darr,points_dir,bif,init)
    
    %set error value to default
    err = 0;
    
    %find and test eigenvalue of first F point
    [e_val(1),e_sol(:,1),e_bif(1),r_flag(1)] = calc_eigenvalue(alpha,F_points(1),L,N,sol);
    if e_val(1) >= 0
        err = 1;
        return
    elseif bif == 2 && e_bif(1) == 1 && r_flag(1) == 1
        err = 1;
        return
    end
    
    %find and test eigenvalue of second F point
    [e_val(2),e_sol(:,2),e_bif(2),r_flag(2)] = calc_eigenvalue(alpha,F_points(2),L,N,e_sol(:,1));
    if e_val(2) >= 0
        err = 2;
        return
    elseif bif == 2 && e_bif(2) == 1 && r_flag(2) == 1
        err = 2;
        return
    end
    
    %find and test eigenvalue of third F point
    [e_val(3),e_sol(:,3),e_bif(3),r_flag(3)] = calc_eigenvalue(alpha,F_points(3),L,N,e_sol(:,2));
    if bif == 1 && e_val(3) <= 0
        err = 3;
        return;
    elseif bif == 2 && ((e_val(3) < e_val(1) && e_val(1)-e_val(3) > 8*abs(e_val(1))) ...
                    ||  (e_val(3) < e_val(2) && e_val(2)-e_val(3) > 5*abs(e_val(2))) || e_val(3) >= 0)
        err = 3;
        return;
    elseif bif == 2 && e_bif(3) == 1 && r_flag(3) == 1
        err = 3;
        return
    end
    
end