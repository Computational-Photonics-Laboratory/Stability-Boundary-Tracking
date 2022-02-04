function [F_calc,F_newpoints,F_sol,e_val,cor_flag] = boundary_locate_dec(alpha,dF,F_points,Fdir,L,N,sol,bif,e_val,darr)

    %if this function is called, largest eigenvalue is moving away from
    %0 as we move towards where we think the boundary is
    %we cant use the polynomial, just have to use a stronger search
    %function to try to find the boundary

    F_calc = -1;
    F_newpoints = F_points;
    F_sol = sol;
    cor_flag = 0;
    
    %starting jump distance
    dis = 10;
    jump_val = dis^4;
    %number of jumps
    num_jumps = 0;

    %start at third point and search towards where the boundary should be
    F_val = F_points(3);
    F_eig = e_val(3);
    ic(:,1) = sol; temp_ic = ic(:,1);
    
    %main loop
    num_loops = 3;
    pos = 1;
    pos_max = 100;
    for m = 1:num_loops

        while F_eig(pos,1) < 0 && pos <= pos_max
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
        if m == num_loops
            break
        else
            jump_val = sqrt(jump_val);
        end
    end
    
    %check to see if final F position is not close to actual boundary
    F_check = F_val(pos) + Fdir * dF;
    [F_eig_check,~,~,~] = calc_eigenvalue(alpha,F_check,L,N,ic(:,pos));
    if F_eig_check < 0 && F_eig_check < F_eig(pos,1) && pos == pos_max
        cor_flag = 1;
        return
    end
    
    %convert to F_calc
    F_calc = F_val(pos);
    
    %recalculate F point set if necessary if boundary is too far away
    if abs(F_calc - F_points(3)) >= dF
        F_newpoints = [(F_calc - Fdir*dF*5/2) (F_calc - Fdir*dF*3/2) (F_calc - Fdir*dF*1/2)];
        [e_val_t,e_sol_t,err_t] = check_eigenvalue(alpha,F_newpoints,L,N,sol,darr,Fdir,bif);
        while err_t ~= 0
            F_newpoints = F_newpoints - Fdir * dF / 2;
            [e_val_t,e_sol_t,err_t] = check_eigenvalue(alpha,F_newpoints,L,N,sol,darr,Fdir,bif);
        end
        [~,F_sol,~,~] = calc_eigenvalue(alpha,F_newpoints(1),L,N,F_sol);
    end
    
end
