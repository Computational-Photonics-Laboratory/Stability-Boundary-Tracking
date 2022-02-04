function [e_val,e_sol,err] = check_eigenvalue_polar(theta_points,theta_dir,darr,L,N,sol,r,ref_alpha,ref_F,bif)
    
    %set error value to default
    err = 0;
    e_val = [0 0 0];
    
    for n = 1:3
        alpha_calc(n) = ref_alpha + r * cos(theta_points(n));
        F_calc(n) = ref_F + r * sin(theta_points(n));
    end
    
    %find and test eigenvalue of first F point
    [e_val(1),e_sol(:,1),e_bif(1),r_flag(1)] = calc_eigenvalue(alpha_calc(1),F_calc(1),L,N,sol);
    if e_val(1) >= 0
        err = 1;
        return
    elseif bif == 2 && e_bif(1) == 1 && r_flag(1) == 1
        err = 1;
        return
    end
    
    %find and test eigenvalue of second F point
    [e_val(2),e_sol(:,2),e_bif(2),r_flag(2)] = calc_eigenvalue(alpha_calc(2),F_calc(2),L,N,e_sol(:,1));
    if e_val(2) >= 0
        err = 2;
        return
    elseif bif == 2 && e_bif(2) == 1 && r_flag(2) == 1
        err = 2;
        return
    end
    
    %find and test eigenvalue of third F point
    [e_val(3),e_sol(:,3),e_bif(3),r_flag(3)] = calc_eigenvalue(alpha_calc(3),F_calc(3),L,N,e_sol(:,2));
    if bif == 1 && e_val(3) <= 0
        err = 3;
        return;
    elseif bif == 2 && e_val(3) >= 0
        err = 3;
        return;
    elseif bif == 2 && e_bif(3) == 1 && r_flag(3) == 1
        err = 3;
        return
    end
    
end