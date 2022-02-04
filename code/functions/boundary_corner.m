function [alpha_bound,F_bound,theta_sol,alpha_diff,Fdir,end_flag,calc_err,bif] = boundary_corner(alpha,da,adir,F,Fdir,alpha_check,F_check,alim_h,alim_l,L,N,sol,slope_tol,dtol)

    disp('corner method called')

    end_flag = 0;
    calc_err = 0;
    
    %number of initial points to be calculated
    num_points = 3;
    
    %determine angular direction to sweep
    theta_dir = adir / Fdir;
    
    %discretization of unit circle
    step_num = 100; step_num_0 = step_num;
    step_inc = step_num / 4;
    step_lim = 3.5*step_num;
    
    %distance from corner to nearest point, corner to farther point
    range_min = sqrt( (alpha(end)-alpha(end-1))^2 + (F(end)-F(end-1))^2 );
    range_max = sqrt( (alpha(end)-alpha(1))^2 + (F(end)-F(1))^2 );
    for n = 1:length(alpha)-1
        r(n) = sqrt( (alpha(end)-alpha(end-n))^2 + (F(end)-F(end-n))^2 );
    end
    
    offset_lim = floor(max(r) / min(r));
    dis = min(r);
    if offset_lim < 2*length(alpha)
        offset_lim = 2*length(alpha);
        dis = (max(r) - min(r)) / offset_lim;
    end
    
    %starting offset value
    if offset_lim < step_lim / step_inc
        offset_max = offset_lim;
    else
        offset_max = round(offset_lim * 0.4);
    end
    %discretization of remaining portion
    offset_inc = ceil((offset_lim-offset_max) * (step_inc / step_lim));
    
    offset = 0;
    offset_flag = 0;
    n = 1;
    flag = 1;
    %create circle of points around corner per nearby point
    alpha_store = []; F_store = []; theta_store = [];
    while flag == 1
        while n <= num_points
            range = min(r) + dis * (n + offset);
            [alpha_next(1:2,:),F_next(1:2,:),theta_next(1:2,:),pos1(n),pos2(n),dtheta(n,1)] = find_nextboundary(alpha,F,range,r,step_num,theta_dir);
            [stab_indices(n,:),bif(n),sol_init(:,n),err] = find_circlestability(alpha_next(1,:),F_next(1,:),L,N,sol(:,end-pos1(n)));
            alpha_store(n,:) = alpha_next(1,:); F_store(n,:) = F_next(1,:); theta_store(n,:) = theta_next(1,:);
            
            if err == 1 && pos2(n) ~= 0
                [stab_indices(n,:),bif(n),sol_init(:,n),err] = find_circlestability(alpha_next(2,:),F_next(2,:),L,N,sol(:,end-pos2(n)));
                alpha_store(n,:) = alpha_next(2,:); F_store(n,:) = F_next(2,:); theta_store(n,:) = theta_next(2,:);
            end
            
            if err == 1
                offset = offset + 1;
            else
                n = n + 1;
            end
            %end of alpha/F corner array exit condition
            if dis == min(r) && n <= num_points
                if n + offset + 1 > offset_max
                    offset_flag = 1;
                    break
                end
            else
                if n + offset > offset_max
                    offset_flag = 1;
                    break
                end
            end
        end
        %if error, change number of steps in unit circle discretization and
        %try again - will add exit condition to this eventually
        if offset_flag == 1
            %increase discretization value
            step_num = step_num + step_inc;
            %increase maximum distance from corner
            if offset_max + offset_inc <= offset_lim
                offset_max = offset_max + offset_inc;
            elseif offset_max + offset_inc > offset_lim && offset_max < offset_lim
                offset_max = offset_lim;
            end
            n = 1;
            offset = 0;
            offset_flag = 0;
            alpha_next = []; F_next = []; theta_next = [];
            alpha_store = []; F_store = []; theta_store = [];
            if step_num > step_lim
                disp("error calculating points after corner")
                %set specific error value to 1, then return unchanged
                %variables
                calc_err = 1;
                alpha_bound = alpha;
                F_bound = F;
                theta_sol = sol;
                alpha_diff = 0;
                bif = 0;
                return
            end
        else
            flag = 0;
        end
    end
    
    %create point sets
    for n = 1:num_points
        alpha_points(n,:) = [alpha_store(n,stab_indices(n,1)) alpha_store(n,stab_indices(n,2)) alpha_store(n,stab_indices(n,3))];
        F_points(n,:) = [F_store(n,stab_indices(n,1)) F_store(n,stab_indices(n,2)) F_store(n,stab_indices(n,3))];
        theta_points(n,:) = [theta_store(n,stab_indices(n,1)) theta_store(n,stab_indices(n,2)) theta_store(n,stab_indices(n,3))];
    end
    
    
    %get first three points after corner
    for n = 1:num_points
        dis(n,1) = sqrt((alpha(end)-alpha_points(n,1))^2 + (F(end)-F_points(n,1))^2);
        [theta_calc(n,1),err(n,1),theta_points(n,:),theta_sol(:,n),e_val(n,:),cor_flag] = boundary_locate_polar(theta_points(n,:),theta_dir,...
            dtheta(n),L,N,sol_init(:,n),dis(n,1),alpha(end),F(end),bif(n));
        
        num_cor = 1; max_cor = 10;
        err_prev = err(n,1);
        %modify theta point and repeat if a point is in the wrong location
        while err(n,1) ~= 0 && num_cor <= max_cor && cor_flag == 0
            %fix point set and run again
            theta_points(n,:) = locate_check_polar(theta_points(n,:),dtheta(n),err(n,1),theta_dir,bif(n));
            [theta_calc(n,1),err(n,1),theta_points(n,:),theta_sol(:,n),e_val(n,:),cor_flag] = boundary_locate_polar(theta_points(n,:),theta_dir,...
                dtheta(n),L,N,theta_sol(:,n),dis(n,1),alpha(end),F(end),bif(n));
            if err(n,1) == err_prev
                num_cor = num_cor + 1;
            elseif err(n,1) ~= err_prev && err(n,1) ~= 0
                err_prev = err(n,1);
                num_cor = 0;
            end
        end
        
        alpha_bound(n) = alpha(end) + dis(n,1) * cos(theta_calc(n,1));
        F_bound(n) = F(end) + dis(n,1) * sin(theta_calc(n,1));
        
        fprintf('alpha = %15.10f calculated\n',alpha_bound(n));
        
        %calculate dF/dalpha slope
        if n == 1
            F_diff = F_bound(n) - F(end);
            alpha_diff = alpha_bound(n) - alpha(end);
            slope(n) = abs(F_diff/alpha_diff);
        else
            F_diff = F_bound(n) - F_bound(n-1);
            alpha_diff = alpha_bound(n) - alpha_bound(n-1);
            slope(n) = abs(F_diff/alpha_diff);
        end
    end

    
    
    %continue with corner method until normal method can be resumed
    d_slope = abs(slope(n)-slope(n-1));    
    %starting distance
    dis_0 = dis(n,1) - dis(n-1,1);
    dis_add = dis_0;
    num_extra = 0; extra_lim = 1;
    while d_slope > dtol || slope(n) > slope_tol || num_extra < extra_lim
        
        %exit condition for a larger slope if slope change is smaller
        %kind of arbitrary
        if d_slope < dtol / 5 && slope(n) < slope_tol * 1.5 && num_extra >= extra_lim
            break
        end
        
        %increment n
        n = n + 1;
        num_extra = num_extra + 1;

        %increase dis if d_slope is small enough
        if d_slope < dtol / 2 && dis_add < dis_0*8
            dis_add = dis_add * 2;
        elseif d_slope > 2*dtol
            dis_add = dis_0;
        end
        dis(n,1) = dis(n-1,1) + dis_add;
        
        %calculate slope from two previous boundary points
        theta_slope = (theta_calc(n-1,1) - theta_calc(n-2,1)) / dis(n-1,1);
        %calculate nexst boundary point set using slope
        theta_points(n,:) = theta_points(n-1,:) + theta_slope * dis(n,1);
        
        dtheta(n,1) = dtheta(n-1,1);
        bif(n) = bif(n-1);
        
        %recalculate next set of points - using same algorithm as previous
        [theta_calc(n,1),err(n,1),theta_points(n,:),theta_sol(:,n),e_val(n,:),cor_flag] = boundary_locate_polar(theta_points(n,:),theta_dir,...
            dtheta(n),L,N,theta_sol(:,n-1),dis(n,1),alpha(end),F(end),bif(n-1));
        num_cor = 1; max_cor = 10;
        err_prev = err(n,1);
        %modify theta point and repeat if a point is in the wrong location
        while err(n,1) ~= 0 && num_cor <= max_cor && cor_flag == 0
            %fix point set and run again
            theta_points(n,:) = locate_check_polar(theta_points(n,:),dtheta(n),err(n,1),theta_dir,bif(n));
            [theta_calc(n,1),err(n,1),theta_points(n,:),theta_sol(:,n),e_val(n,:),cor_flag] = boundary_locate_polar(theta_points(n,:),theta_dir,...
                dtheta(n),L,N,theta_sol(:,n),dis(n,1),alpha(end),F(end),bif(n-1));
            if err(n,1) == err_prev
                num_cor = num_cor + 1;
            elseif err(n,1) ~= err_prev && err(n,1) ~= 0
                err_prev = err(n,1);
                num_cor = 0;
            end
        end
        if theta_calc(n,1) == -1
            break
        end
        
        alpha_bound(n) = alpha(end) + dis(n,1) * cos(theta_calc(n,1));
        F_bound(n) = F(end) + dis(n,1) * sin(theta_calc(n,1));
        
        fprintf('alpha = %15.10f calculated\n',alpha_bound(n));
        
        %recalculate slope and change in slope
        F_diff = F_bound(n) - F_bound(n-1);
        alpha_diff = alpha_bound(n) - alpha_bound(n-1);
        slope(n) = abs(F_diff/alpha_diff);
        d_slope = abs(slope(n)-slope(n-1));
        
        %check for reaching starting boundary point exit condition
        %check if latest point matches exit boundary point
        %check alpha first - if alpha_check is between last step values
        if (alpha_diff > 0 && alpha_bound(n) > alpha_check && alpha_bound(n-1) < alpha_check) || ...
           (alpha_diff < 0 && alpha_bound(n) < alpha_check && alpha_bound(n-1) > alpha_check)
            %check F now - if F_check is between last step values
            if (F_diff > 0 && F_bound(n) > F_check && F_bound(n-1) < F_check) || ...
               (F_diff < 0 && F_bound(n) < F_check && F_bound(n-1) > F_check)
                alpha_bound(n) = alpha_check;
                F_bound(n) = F_check;
                Fdir = 0;
                bif = 0;
                end_flag = 1;
                disp("reached starting boundary_point.");
                return
            end
        end
        
        if alpha_bound(n) > alim_h
            end_flag = 2;
            disp("reached upper alpha limit")
            return
        elseif alpha_bound(n) < alim_l
            end_flag = 2;
            disp("reached lower alpha limit")
            return
        end
        
    end
    
    %recalculate Fdir for next smooth section
    for m = 1:3
        F_points_t(m) = F(end) + dis(n,1) * sin(theta_points(n,m));
    end
    if F_points_t(3) - F_points_t(2) > 0
        Fdir = 1;
    else
        Fdir = -1;
    end
    
    %return bifurcation type
    bif = bif(n);
    
    disp('corner method exited')

end