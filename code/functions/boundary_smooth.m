function [alpha_bound,F_bound,sol_bound,num_corner,F_points_corner,end_flag] = boundary_smooth(alpha,da,adir,F_points,dF,Fdir,F_bound_curr,alpha_check,F_check,alim_h,alim_l,L,N,sol,bif,dtol,slope_init,grow_lim)

    end_flag = 0;

    alpha_bound = []; F_bound = []; sol_bound = [];
    
    %intial alpha point
    alpha_0 = alpha;
    
    alpha_lim = alpha;

    %dalpha and dF values
    da_0 = da; dF_0 = dF;
    da_rat = da / da_0;
    num_divis = 0; 
    divis_limit = grow_lim + 2;
    F_steps_max = 2*(divis_limit);

    %loop variables
    %number of boundary points calculated
    num_points = 0;
    %stability flag
    stable = 1;
    %counter to see if decreased recently - will not increase next chance
    %to if on
    dec_flag = 0;
    %variables to store point to be started from if recalculating the
    %corner
    num_corner = 0;
    F_points_corner = [0; 0; 0];
    
    %initial F next point set
    F_points_next = F_points;
    F_points_prev = F_points;


    %main loop start - repeats until function ends
    while 1
        
        %save previous solution
        prev_sol = sol;
        
        %increment boundary point counter
        num_points = num_points + 1;
    
        while num_divis < divis_limit

            %increment alpha
            alpha = round(alpha + adir*da,8);
            
            %stop alpha from increasing further if da is below the limit
            if da < da_0 / 2^grow_lim && round(alpha,8) == round(alpha_lim,8)
                da = da / 2;
                da_rat = da / da_0;
                dec_flag = 1;
                alpha = alpha - adir*da;
                num_divis = num_divis + 1;
                F_points = F_points_prev + slope(num_points-1) * da;
            end

    
            %calculate boundary point at new alpha
            %changed to check in appropriate alpha direction
            %use solution of last boundary point calculated at initial condition
            if num_points == 1
                [F_bound_next,err,F_points,sol,e_val,cor_flag] = boundary_locate(da_rat,alpha,dF,F_points,Fdir,L,N,sol,bif,0);
            else
                [F_bound_next,err,F_points,sol,e_val,cor_flag] = boundary_locate(da_rat,alpha,dF,F_points,Fdir,L,N,sol,bif,0);
            end
    
            %check if error while calculating boundary
            num_cor = 1;
            err_prev = err;
            while num_cor <= F_steps_max-2*num_divis && err ~= 0 && cor_flag == 0
                if bif == 1
                    e_val_prev = e_val(err);
                end
                %move point set 
                [F_points] = locate_check(dF,Fdir,err,F_points,bif);
                %call function again
                if num_points == 1
                    [F_bound_next,err,F_points,sol,e_val,cor_flag] = boundary_locate(da_rat,alpha,dF,F_points,Fdir,L,N,sol,bif,0);
                else
                    [F_bound_next,err,F_points,sol,e_val,cor_flag] = boundary_locate(da_rat,alpha,dF,F_points,Fdir,L,N,sol,bif,0);
                end
                %check error results if any
                if err == err_prev && cor_flag == 0
                    num_cor = num_cor + 1;
                    if bif == 1
                        if err ~= 3 && e_val_prev > 0 && e_val(err) > e_val_prev
                            break
                        elseif err == 3 && e_val_prev < 0 && e_val(err) < e_val_prev
                            break
                        end
                    end 
                elseif err ~= err_prev && err ~= 0
                    err_prev = err;
                    num_cor = 0;
                end
            end
            %exit condition - if F_bound_next is found successfully, exit
            %if not, ie = -1, then move back alpha, reduce da, and try again
            if F_bound_next ~= -1
                break
            else
                if num_points == 1
                    %save alpha value failed at
                    alpha_lim = alpha;
                    %reset alpha to previous value
                    alpha = alpha - adir*da;
                    da = da / 2;
                    da_rat = da / da_0;
                    dec_flag = 1;
                    num_divis = num_divis + 1;
                    %recalculate F_points based on different da
                    F_points = F_points - slope_init * da;
                else
                    %save alpha value failed at
                    alpha_lim = alpha;
                    %reset alpha to previous value
                    alpha = alpha - adir*da;
                    %change da
                    if da > da_0
                        da = da_0;
                        da_rat = 1;
                        dec_flag = 1;
                    else
                        da = da / 2;
                        da_rat = da / da_0;
                        dec_flag = 1;
                        num_divis = num_divis + 1;
                        %save values before decreasing for corner
                        %calculation if first decrease in da
                        if num_divis == 1
                            if num_points == 2
                                num_corner = num_points - 1;
                            else
                                num_corner = num_points - 2;
                            end
                            F_points_corner = F_points_bound(num_corner,:) + slope(num_corner) * da;
                        end
                    end
                    %recalculate F_points based on different da  
                    F_points = F_points_prev + slope(num_points-1) * da;               
                end
            end
    
        end
        %exit condition if failed to find a good boundary point
        if F_bound_next == -1
            break
        end

    
        %calculate slope from curr and next boundary points
        slope(num_points) = (F_bound_next - F_bound_curr) / da;
    
        %calculate change in slope
        if num_points == 1
            d_slope = abs(slope(num_points));
        else
            d_slope = abs(slope(num_points)-slope(num_points-1));
        end
        
        %decrement da if slope change is too hight. if just decremented da, 
        %dont increase da, otherwise increase da if slope change is small
        if dec_flag == 1
            dec_flag = 0;
        elseif da >= da_0/2^grow_lim && da < 4*da_0 && d_slope < dtol/2
            if da < da_0
                num_divis = num_divis - 1;
            end
            da = da * 2;
            da_rat = da / da_0;
        end
        
        
        fprintf('alpha = %15.10f calculated\n',alpha);
        
        %change vars for next iteration
        F_points_next = F_points + slope(num_points) * da;
        %store alpha and F boundary point values
        alpha_bound(num_points) = alpha;
        F_bound(num_points) = F_bound_next;
        F_points_bound(num_points,:) = F_points;
        sol_bound(:,num_points) = sol;
        %save previous point set
        F_points_prev = F_points;
        %change F_points to be last two calculated point set
        F_points = F_points_next;
        %change F_bound_curr to be F_bound_next
        F_bound_curr = F_bound_next;
        
        %number of divisions of da exit condition
        if num_divis >= divis_limit
            break
        end
        
        %check if latest point matches exit boundary point
        %check alpha first - if alpha_check is inside range of next da step
        if (adir == 1 && alpha < alpha_check && alpha + da > alpha_check) || (adir == -1 && alpha > alpha_check && alpha - da < alpha_check)
            %check F now - 3*dF minimum range
            if da > da_0
                da_rat = da / da_0;
            else
                da_rat = 1;
            end
            F_proj = F_bound_next + slope(num_points) * abs(alpha_check - alpha);
            if abs(F_proj - F_check) < (3*dF) * da_rat
                alpha_bound = [alpha_bound alpha_check];
                F_bound = [F_bound F_check];
                end_flag = 1;
                disp("reached starting boundary_point.");
                return
            end
        end
        
        if alpha > alim_h
            end_flag = 2;
            disp("reached upper alpha limit")
            return
        elseif alpha < alim_l
            end_flag = 2;
            disp("reached lower alpha limit")
            return
        end
    
    end

end