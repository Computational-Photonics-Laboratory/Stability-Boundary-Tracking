function [F_calc,err,F_newpoints,F_sol,e_val,cor_flag] = boundary_locate(da_rat,alpha,dF,F_points,Fdir,L,N,sol,bif_prev,init)

    F_calc = -1;
    err = 0;
    F_newpoints = F_points;
    F_sol = sol;
    e_val = [0 0 0];
    cor_flag = 0;
    bif = bif_prev;
    
    arr_scale = 1e-3; darr = dF*arr_scale;
    %caculate largest real eigenvalue parts for F_point set, return an
    %error value if necessary depending on the point that caused the error
    [e_val,e_sol,err] = check_eigenvalue(alpha,F_points,L,N,sol,darr,Fdir,bif,init);
    if err ~= 0
        return
    end
    
    thresh = 1e-6; %should be same as threshold used in calc_eigenvalue for largest real eigenvalue
    if e_val(1) > e_val(2) && e_val(2) > e_val(3) && e_val(2) - e_val(3) > 5*thresh
        [F_calc,F_newpoints,F_sol,e_val,cor_flag] = boundary_locate_dec(alpha,dF,F_points,Fdir,L,N,e_sol(:,3),bif,e_val,darr);
        if cor_flag == 1
            F_sol = sol;
        end
        return
    end
    
    %set solution to return as solution to innermost F point
    F_sol = e_sol(:,1);
    
    poly = polyfit(F_points,e_val,2);
    %create polynomial
    if bif == 1
        if Fdir == 1
            arr = F_points(1):darr:F_points(3);
        else
            arr = F_points(3):darr:F_points(1);
        end
    elseif bif == 2
        if Fdir == 1
            arr = F_points(1):darr:F_points(3)+dF;
        else
            arr = F_points(3)-dF:darr:F_points(1);
        end       
    end
    res = polyval(poly,arr);
    %find and store F value with corresponding e value closest to 0
    e_min = 1;
    count = 1;
    for n = 1:length(res)
        %check if closer to 0 then previous stored smallest value
        if abs(res(n)) < abs(e_min) && res(n) < 0
            count = n;
            e_min = res(n);
            F_val = arr(n);
        end
    end
    
    
    if bif == 1
    
    %current jump distance
    jump_val = ceil(sqrt(1/arr_scale));
    %number of jumps
    num_jumps = 0;

    %get initial eigenvalue
    [F_eig,ic(:,1),~,~] = calc_eigenvalue(alpha,F_val,L,N,e_sol(:,2));
    
    %if greater than 0, we want to rapidly move it back until less than 0
    %moving F away from boundary location
    if F_eig >= 0
        %initial condition
        ic(:,1) = e_sol(:,2);
        while F_eig >= 0
            %jump
            F_val = F_val - Fdir * (darr * jump_val);
            %check if F_val has gone out of bounds
            if F_val * Fdir < F_points(2) * Fdir
                F_val = F_points(2);
            end
            num_jumps = num_jumps + 1;
            %recalculate F_eig
            [F_eig,temp_ic,~,~] = calc_eigenvalue(alpha,F_val,L,N,ic(:,1));
        end
        %update initial condition for new starting F_val
        ic(:,1) = temp_ic;
    else
        %F_eig already negative - no movement necessary
        %set temporary initial condition variable
        temp_ic = ic(:,1);
    end
    
    %main loop
    num_loops = 3;
    pos = 1;
    for m = 1:num_loops
        %jump loop
        while F_eig(pos,1) < 0
            num_jumps = num_jumps + 1;
            pos = pos + 1;
            %update initial condition
            ic(:,pos-1) = temp_ic;
            %recalculate F_eig
            F_val(pos,1) = F_val(pos-1,1) + Fdir * (darr * jump_val);
            [F_eig(pos,1),temp_ic,~,~] = calc_eigenvalue(alpha,F_val(pos,1),L,N,ic(:,pos-1));
        end
        %change back variables to last good F_val result
        pos = pos - 1;
        temp_ic = ic(:,pos);
        %decrease jump val - set to 1 if on last loop
        if m >= num_loops - 1
            jump_val = 1;
        elseif jump_val == 1
            break;
        else
            jump_val = ceil(sqrt(jump_val));
        end
    end
    
    %convert to F_calc
    F_calc = F_val(pos);
    
    end
    
    
    
    if bif == 2
    
    %starting jump distance
    dis = 8;
    jump_val = dis^4;
    %number of jumps
    num_jumps = 0;

    %get initial eigenvalue, set initial condition and search ,direction
    [F_eig,ic(:,1),~,~] = calc_eigenvalue(alpha,F_val,L,N,e_sol(:,3));
    
    %if greater than 0, we want to rapidly move it back until less than 0
    %moving F away from boundary location
    if F_eig >= 0
        diff = (F_val - F_points(3)) / 10;
        while F_eig >= 0
            %jump
            if F_val * Fdir > F_points(3) * Fdir
                F_val = F_val - diff;
            else
                F_val = F_points(3);
                F_eig = e_val(3);
                temp_ic = e_sol(:,3);
                break
            end
            num_jumps = num_jumps + 1;
            %recalculate F_eig
            [F_eig,temp_ic,~,~] = calc_eigenvalue(alpha,F_val,L,N,e_sol(:,3));
        end
        %update initial condition for new starting F_val
        ic(:,1) = temp_ic;
    else
        %F_eig already negative - no movement necessary
        %set temporary initial condition variable
        temp_ic = ic(:,1);
    end
    
    
    %check adjacent F_eig values to determine if on wrong side of maximum
    %real eigenvalue part
    test_jump = dis*2;
    recalc_flag = 0;
    F_test = F_val - Fdir * test_jump * darr;
    [F_eig_n,~,~,~] = calc_eigenvalue(alpha,F_test,L,N,ic(:,1));
    %track away from the boundary location until on correct side
    while F_eig_n > F_eig || F_eig >= 0 || F_eig_n >= 0
        if F_val * Fdir < F_points(1) * Fdir
            cor_flag = 1;
            F_sol = sol;
            return
        end
        recalc_flag = 1;
        F_val = F_test;
        F_test = F_val - Fdir * test_jump * darr;
        F_eig = F_eig_n;
        [F_eig_n,~,~,~] = calc_eigenvalue(alpha,F_test,L,N,ic(:,1));
    end
    
    
    %main loop
    num_loops = 3;
    pos = 1;
    pos_max = 100;
    for m = 1:num_loops
        %jump loop
        eig_diff = 1;
        while F_eig(pos,1) < 0 && eig_diff > 0 && pos <= pos_max
            num_jumps = num_jumps + 1;
            pos = pos + 1;
            %update initial condition
            ic(:,pos-1) = temp_ic;
            %recalculate F_eig
            F_val(pos,1) = F_val(pos-1,1) + Fdir * (darr * jump_val);
            [F_eig(pos,1),temp_ic,~,~] = calc_eigenvalue(alpha,F_val(pos,1),L,N,ic(:,pos-1));
            %check change in largest real eigenvalue part
            eig_diff = F_eig(pos,1) - F_eig(pos-1,1);
        end
        %change back variables to last good F_val result
        pos = pos - 1;
        temp_ic = ic(:,pos);
        if m == num_loops
            break
        else
            jump_val = sqrt(jump_val);
        end
        
        %recheck F direction value
        if pos > 1 && pos < pos_max
            [F_eig_n,~,~,~] = calc_eigenvalue(alpha,F_val(pos,1)-Fdir*(darr*test_jump),L,N,ic(:,pos-1));
            if F_eig(pos,1) - F_eig_n <= 0
                pos = pos - 1;
                temp_ic = ic(:,pos);
            end
        end
    end
    
    %check to see if final F position is not close to actual boundary
    F_check = F_val(pos) + Fdir * dF;
    [F_eig_check,~,~,~] = calc_eigenvalue(alpha,F_check,L,N,ic(:,pos));
    if F_eig_check < 0 && F_eig_check > F_eig(pos,1) && pos == pos_max
        cor_flag = 1;
        F_sol = sol;
        return
    end
    
    %convert to F_calc
    F_calc = F_val(pos);
    
    %recalculate F point set if necessary - if boundary is too far  away or
    %inside previous F point set
    if abs(F_calc - F_points(3)) >= dF || recalc_flag == 1
        F_newpoints = [(F_calc - Fdir*dF*5/2) (F_calc - Fdir*dF*3/2) (F_calc - Fdir*dF*1/2)];
        [~,~,err_t] = check_eigenvalue(alpha,F_newpoints,L,N,sol,darr,Fdir,bif);
        while err_t ~= 0
            F_newpoints = F_newpoints - Fdir * dF / 2;
            [~,~,err_t] = check_eigenvalue(alpha,F_newpoints,L,N,sol,darr,Fdir,bif);
        end
        [~,F_sol,~,~] = calc_eigenvalue(alpha,F_newpoints(1),L,N,F_sol);   
    end
    
    end
    
end
