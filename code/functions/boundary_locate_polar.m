function [theta_calc,err,theta_newpoints,theta_sol,e_val,cor_flag] = boundary_locate_polar(theta_points,theta_dir,dtheta,L,N,sol,r,ref_alpha,ref_F,bif_prev)

    %bif determines behavior of function - hopf or saddle node

    theta_calc = -1;
    err = 0;
    theta_newpoints = theta_points;
    theta_sol = sol;
    e_val = [0 0 0];
    cor_flag = 0;
    bif = bif_prev;
    
    arr_scale = 1e-3; darr = dtheta*arr_scale;
    %caculate largest real eigenvalue parts for theta_point set, return an
    %error value if necessary depending on the point that caused the error
    [e_val,e_sol,err] = check_eigenvalue_polar(theta_points,theta_dir,darr,L,N,sol,r,ref_alpha,ref_F,bif);
    if err ~= 0
        return
    end
    
    %set solution to return as solution to innermost theta point
    theta_sol = e_sol(:,1);
    poly = polyfit(theta_points,e_val,2);
    %create polynomial
    if bif == 1
        if theta_dir == 1
            arr = theta_points(1):darr:theta_points(3);
        else
            arr = theta_points(3):darr:theta_points(1);
        end
    elseif bif == 2
        if theta_dir == 1
            arr = theta_points(1):darr:theta_points(3)+dtheta;
        else
            arr = theta_points(3)-dtheta:darr:theta_points(1);
        end      
    end
    res = polyval(poly,arr);
    %find and store theta value with corresponding e value closest to 0
    e_min = 1;
    count = 1;
    for n = 1:length(res)
        if abs(res(n)) < abs(e_min) && res(n) < 0
            count = n;
            e_min = res(n);
            theta_val = arr(n);
        end
    end
    
    
    if bif == 1
    
    len = 1/arr_scale;
    %current jump distance
    jump_val = ceil(sqrt(len));
    %number of jumps
    num_jumps = 0;

    %get initial eigenvalue
    alpha_calc = ref_alpha + r * cos(theta_val);
    F_calc = ref_F + r * sin(theta_val);
    [theta_eig,ic(:,1),~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,e_sol(:,2));
    
    %if greater than 0, we want to rapidly move it back until less than 0
    %moving theta away from boundary location
    if theta_eig >= 0
        %initial condition
        ic(:,1) = e_sol(:,2);
        while theta_eig >= 0
            %jump
            theta_val = theta_val - theta_dir * (darr * jump_val);
            %check if theta_val has gone out of bounds
            if theta_val * theta_dir < theta_points(2) * theta_dir
                theta_val = theta_points(2);
            end
            num_jumps = num_jumps + 1;
            %recalculate theta_eig
            alpha_calc = ref_alpha + r * cos(theta_val);
            F_calc = ref_F + r * sin(theta_val);
            [theta_eig,temp_ic,~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,ic(:,1));
        end
        %update initial condition for new starting theta_val
        ic(:,1) = temp_ic;
    else
        %theta_eig already negative - no movement necessary
        %set temporary initial condition variable
        temp_ic = ic(:,1);
    end
    
    %main loop
    num_loops = 3;
    pos = 1;
    for m = 1:num_loops
        %jump loop
        while theta_eig(pos,1) < 0
            num_jumps = num_jumps + 1;
            pos = pos + 1;
            %update initial condition
            ic(:,pos-1) = temp_ic;
            %recalculate theta_eig
            theta_val(pos,1) = theta_val(pos-1,1) + theta_dir * (darr * jump_val);
            alpha_calc = ref_alpha + r * cos(theta_val(pos,1));
            F_calc = ref_F + r * sin(theta_val(pos,1));
            [theta_eig(pos,1),temp_ic,~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,ic(:,pos-1));
        end
        %change back variables to last good theta_val result
        pos = pos - 1;
        temp_ic = ic(:,pos);
        %decrease jump val - set to 1 if on last loop
        if m == num_loops
            break
        elseif m == num_loops - 1
            jump_val = 1;
        else
            jump_val = ceil(sqrt(jump_val));
        end
    end
    
    %convert to theta_calc
    theta_calc = theta_val(pos);
    
    end
    
    
    
    if bif == 2
    
    %create new array of that represents the searchable area
    len = 1/arr_scale;
    %starting jump distance
    dis = 8;
    jump_val = dis^4;
    %number of jumps
    num_jumps = 0;
    
    %get initial eigenvalue, set initial condition and search ,direction
    alpha_calc = ref_alpha + r * cos(theta_val);
    F_calc = ref_F + r * sin(theta_val);
    [theta_eig,ic(:,1),~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,e_sol(:,3));
    
    %if theta_eig greater than 0, rapidly move theta_val towards point set until less than 0
    if theta_eig >= 0
        dis = (theta_val - theta_points(3)) / 10;
        while theta_eig >= 0
            %jump
            if theta_val * theta_dir > theta_points(3) * theta_dir
                theta_val = theta_val - dis;
            else
                theta_val = theta_points(3);
                theta_eig = e_val(3);
                temp_ic = e_sol(:,3);
                break
            end
            num_jumps = num_jumps + 1;
            %recalculate theta_eig
            alpha_calc = ref_alpha + r * cos(theta_val);
            F_calc = ref_F + r * sin(theta_val);
            [theta_eig,temp_ic,~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,ic(:,1));
        end
        %update initial condition for new starting theta_val
        ic(:,1) = temp_ic;
    else
        %theta_eig already negative - no movement necessary
        %set temporary initial condition variable
        temp_ic = ic(:,1);
    end
    
    %main loop
    num_loops = 3;
    pos = 1;
    pos_max = 100;
    for m = 1:num_loops
        %jump loop
        eig_diff = 1;
        while theta_eig(pos,1) < 0 && pos <= pos_max
            num_jumps = num_jumps + 1;
            pos = pos + 1;
            %update initial condition
            ic(:,pos-1) = temp_ic;
            %recalculate theta_eig
            theta_val(pos,1) = theta_val(pos-1,1) + theta_dir * (darr * jump_val);
            alpha_calc = ref_alpha + r * cos(theta_val(pos,1));
            F_calc = ref_F + r* sin(theta_val(pos,1));
            [theta_eig(pos,1),temp_ic,~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,ic(:,pos-1));         
        end
        %change back variables to last good theta_val result
        pos = pos - 1;
        temp_ic = ic(:,pos);
        if m == num_loops
            break
        else
            jump_val = sqrt(jump_val);
        end
    end
    
    %convert to theta_calc
    theta_calc = theta_val(pos);
    
    end

    
    if bif == 2
        if abs((theta_calc) - theta_points(3)) >= dtheta
            theta_newpoints = [(theta_calc - theta_dir*dtheta*5/2) (theta_calc - theta_dir*dtheta*3/2) (theta_calc - theta_dir*dtheta*1/2)];
            alpha_calc = ref_alpha + r * cos(theta_newpoints(1));
            F_calc = ref_F + r * sin(theta_newpoints(1));
            [~,~,err_t] = check_eigenvalue_polar(theta_points,theta_dir,darr,L,N,sol,r,ref_alpha,ref_F,bif);
            while err_t ~= 0
                theta_newpoints = theta_newpoints - theta_dir * dtheta / 2;
                alpha_calc = ref_alpha + r * cos(theta_newpoints(1));
                F_calc = ref_F + r * sin(theta_newpoints(1));
                [~,~,err_t] = check_eigenvalue_polar(theta_points,theta_dir,darr,L,N,sol,r,ref_alpha,ref_F,bif);
            end
            [~,theta_sol,~,~] = calc_eigenvalue(alpha_calc,F_calc,L,N,theta_sol);
        end
    end
    
end