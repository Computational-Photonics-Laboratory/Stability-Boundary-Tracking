function [alpha_bound,F_bound,sol,end_flag] = boundary_wrapper(alpha,da,adir,F,dF,Fdir,L,N,uout,alim_h,alim_l,iter)

    sol = uout;
    da_0 = da; dF_0 = dF;

    %largest slope acceptable to the problem
    slope_tol = 2.5;
    %largest slope change that is acceptable to the program
    dtol = 0.6;
    %number of times da_0 can be divided by 2^x before starting to call the 
    %corner method. Can be increased to give more resolution, but will
    %significantly slow down the program
    grow_lim = 3;

    %bifurcation type variable
    bif = 0;

    %call set up function - get two points on boundary in vectors, F point set
    %surrounding the next projected point, and the solution for the last point
    [alpha_loc,F_loc,F_points,sol_loc,bif] = boundary_init(alpha,da,adir,F,dF,Fdir,L,N,uout,iter);

    %alpha and F container vectors - used to hold computed boundary points
    alpha_bound = alpha_loc;
    F_bound = F_loc;
    sol_bound = sol_loc;

    %initial slope of the boundary
    F_slope = (F_loc(end) - F_loc(end-1)) / da;

    %starting boundary point value - used to check if completed boundary
    alpha_check = alpha_bound(1);
    F_check = F_bound(1);

    %iteration counter
    z = 0;

    %main loop
    while 1
    
        %temporary values used for corner recalculation if necessary
        alpha_temp = alpha_bound;
        F_temp = F_bound;
        sol_temp = sol_bound;
        F_points_temp = F_points;
    
        %call smooth segment function - alpha and F values passed in are values
        %to make sure not repeated
        [alpha_loc,F_loc,sol_loc,num_corner,F_points_corner,end_flag] = boundary_smooth(alpha_bound(end),da,adir,F_points,dF,Fdir,F_bound(end),alpha_check,F_check,alim_h,alim_l,L,N,sol_bound(:,end),bif,dtol,F_slope,grow_lim);
        %update alpha and F container vectors
        alpha_bound = [alpha_bound alpha_loc];
        F_bound = [F_bound F_loc];
        sol_bound = [sol_bound sol_loc];
    
        %check if starting point or alpha limit is reached
        if end_flag == 1
            return
        elseif end_flag == 2
            sol = sol_bound(:,1);
            return
        end
    
        %check initial boundary smooth result - needs to generate at least 2
        %points to have a total of 4 in bound arrays. If not, throw error and
        %pick a better starting point
        min_points = 6; 
        if z == 0 && length(alpha_bound) < min_points
            error("Not enough initial boundary points could be generated from initial conditions. Please pick a better starting point.");
        end
    
        %need matrix of length 4 to be passed into boundary corner - if loc
        %arrays not long enough, use last few calculated points
        len = length(alpha_loc);
        max_points = 15;
        if len < min_points
            alpha_loc = alpha_bound(end-min_points+1:end);
            F_loc = F_bound(end-min_points+1:end);
            sol_loc = sol_bound(:,end-min_points+1:end);
        elseif len > max_points
            alpha_loc = alpha_bound(end-max_points+1:end);
            F_loc = F_bound(end-max_points+1:end);
            sol_loc = sol_bound(:,end-max_points+1:end);
        end
        
        %increase grow lim if necessary - usually stuck in a corner and
        %need higher resolution
        if len < 5
            grow_lim = 4;
        else
            grow_lim = 3;
        end
    
        %when done, at a corner point or at starting point
        %call corner function to handle
        [alpha_loc,F_loc,sol_loc,alpha_diff,Fdir,end_flag,err,bif_t] = boundary_corner(alpha_loc,da,adir,F_loc,Fdir,alpha_check,F_check,alim_h,alim_l,L,N,sol_loc,slope_tol,dtol);
    
        %recalculate corner if err ~= 0
        if err ~= 0
            disp("recalculating corner")
            %decrease da for more accurate corner calculation
            da = da/2;
            %move alpha/F vectors back 10 spaces if long smooth curve calcuation
            %reset to before previous smooth curve calculation if otherwise
            if len > max_points
                diff = len - num_corner;
                alpha_bound = alpha_bound(1:end-diff);
                F_bound = F_bound(1:end-diff);
                sol_bound = sol_bound(:,1:end-diff);
                F_points = F_points_corner;
                %case for too many error calculations
                if da < da_0 / 4
                    da = da_0 * 4/5;
                    alpha_bound = alpha_temp;
                    F_bound = F_temp;
                    sol_bound = sol_temp;
                    F_points = F_points_temp;
                end
            else
                if isequal(alpha_bound,alpha_temp)
                    alpha_bound = alpha_bound(1:end-min_points+1);
                    F_bound = F_bound(1:end-min_points+1);
                    sol_bound = sol_bound(:,1:end-min_points+1);
                    if bif == 1
                        F_points = [F_bound(end)-Fdir*dF*3/2 F_bound(end)-Fdir*dF*1/2 F_bound(end)+Fdir*dF*1/2];
                    else
                        F_points = [F_bound(end)-Fdir*dF*5/2 F_bound(end)-Fdir*dF*3/2 F_bound(end)-Fdir*dF*1/2];
                    end
                else
                    alpha_bound = alpha_temp;
                    F_bound = F_temp;
                    sol_bound = sol_temp;
                    F_points = F_points_temp;
                end
                if da < da_0 / 8
                    disp("corner could not be recalculated to allow for accurate modelling with previous data. Please pick a different initial condition and try again.")
                    disp("data recorded in alpha_bound/F_bound container vectors is most likely inaccurate after last corner.")
                    return
                end
            end
        else
        
            %update alpha and F container vectors
            alpha_bound = [alpha_bound alpha_loc];
            F_bound = [F_bound F_loc];
            sol_bound = [sol_bound sol_loc];
        
            da = da_0;
            %check if starting point is reached or alpha limit is reached
            if end_flag == 1
                return
            elseif end_flag == 2
                sol = sol_bound(:,1);
                return
            end
    
            %reset alpha and F direction
            if alpha_loc(end) - alpha_loc(end-1) < 0
                adir = -1;
            else
                adir = 1;
            end
    
            bif = bif_t;
    
            %recalculate F_points
            while 1
                if bif == 1
                    F_points = [F_loc(end)-Fdir*dF*3/2 F_loc(end)-Fdir*dF*1/2 F_loc(end)+Fdir*dF*1/2];
                    [E,temp_sol,bif_test,~] = calc_eigenvalue(alpha_loc(end),F_points(2),L,N,sol_bound(:,end));
                else
                    F_points = [F_loc(end)-Fdir*dF*5/2 F_loc(end)-Fdir*dF*3/2 F_loc(end)-Fdir*dF*1/2];
                    [E,temp_sol,bif_test,~] = calc_eigenvalue(alpha_loc(end),F_points(3),L,N,sol_bound(:,end));
                end
                %check if first F_point is in stable region
                if E > 0
                    Fdir = -1 * Fdir;
                    if bif == 1
                        F_points = [F_loc(end)-Fdir*dF*3/2 F_loc(end)-Fdir*dF*1/2 F_loc(end)+Fdir*dF*1/2];
                        [E,temp_sol,bif_test,~] = calc_eigenvalue(alpha_loc(end),F_points(2),L,N,sol_bound(:,end));
                    else
                        F_points = [F_loc(end)-Fdir*dF*5/2 F_loc(end)-Fdir*dF*3/2 F_loc(end)-Fdir*dF*1/2];
                        [E,temp_sol,bif_test,~] = calc_eigenvalue(alpha_loc(end),F_points(3),L,N,sol_bound(:,end));
                    end
                end
                %exit condition
                if E < 0
                    break
                else
                    %set Fdir back to original value
                    Fdir = -1 * Fdir;
                
                    alpha_loc = alpha_bound(end-min_points-3+1:end);
                    F_loc = F_bound(end-min_points-3+1:end);
                    sol_loc = sol_bound(:,end-min_points-3+1:end);
                
                    %if cant exit, call boundary corner and try again
                    [alpha_loc,F_loc,sol_loc,alpha_diff,Fdir,end_flag,err,bif_t] = boundary_corner(alpha_loc,da,adir,F_loc,Fdir,alpha_check,F_check,alim_h,alim_l,L,N,sol_loc,slope_tol,dtol);
                
                    if ~isequal(alpha_bound(end-min_points-3+1:end),alpha_loc)
                        alpha_bound = [alpha_bound alpha_loc];
                        F_bound = [F_bound F_loc];
                        sol_bound = [sol_bound sol_loc];
                    end
                    da = da_0;
                    %check if starting point is reached or alpha limit is reached
                    if end_flag == 1
                        return
                    elseif end_flag == 2
                        sol = sol_bound(:,1);
                        return
                    end
                    %reset alpha and F direction
                    if alpha_loc(end) - alpha_loc(end-1) < 0
                        adir = -1;
                    else
                        adir = 1;
                    end
                    bif = bif_t;
                end
            
            end
        
            if bif_test == bif
                sol_bound(:,end) = temp_sol;
            end
    
            %calculate slope from curr and next boundary points
            F_slope = (F_loc(end) - F_loc(end-1)) / abs(alpha_diff);
            %calculate alpha + 2*da boundary point set using slope
            F_points = F_points + F_slope * da;
    
        end

    end    %end while loop
    
    
end    %end function


