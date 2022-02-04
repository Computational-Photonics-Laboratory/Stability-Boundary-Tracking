function [alpha_bound,F_bound,F_points,sol,bif] = boundary_init(alpha,da,adir,F,dF,Fdir,L,N,sol,iter)
    
    %loop set up
    
    %if second time calling boundary wrapper, move initial F value 50*dF
    %back to allow for large changes in slope
    if iter == 2
        F = F - Fdir*50*dF;
    end
    
    %save initial variables
    alpha_0 = alpha;
    F_0 = F;
    
    dF_0 = dF;
    %large dF for more rapid movement initially
    s = 10;
    dF = s*dF;

    %flag to keep running or exit loop
    stab_flag = 1;
    
    F_arr = [];
    %number of F values stored per loop
    num = [1;1];
    %temporary solution vector
    uout = zeros(2*N,1);
    %difference in largest real eigenvalue part between loops
    eig_diff = [1;1];

    %two loops, one for each alpha
    for k = 1:2
        %inner loop that changes F
        while stab_flag == 1     
            %call calc_eigenvalue on current F/alpha values
            if num(k) == 1
                [E(num(k),k),uout(:,num(k)),bif(num(k),k),~] = calc_eigenvalue(alpha,F,L,N,sol);
            else
                [E(num(k),k),uout(:,num(k)),bif(num(k),k),~] = calc_eigenvalue(alpha,F,L,N,uout(:,num(k)-1));
            end
            
            %recalculate eig_diff
            if num(k) > 1
                eig_diff(num(k),k) = E(num(k),k) - E(num(k)-1,k);
                if eig_diff(num(k),k) * eig_diff(num(k)-1,k) < 0 && abs(eig_diff(num(k),k)) > abs(E(num(k)-1,k)) * 3
                    stab_flag = 0;
                end
                if num(k) > 2 && eig_diff(num(k),k) * eig_diff(num(k)-2,k) < 0 && abs(E(num(k),k) - E(num(k)-2,k)) > abs(E(num(k)-2,k)) * 5
                    stab_flag = 0;
                end 
            end
            
            %set stability flag to 0 if unstable eigenvalue found or if too
            %big a difference in eig_diff and it changes direction
            if E(num(k),k) > 0
                stab_flag = 0;
            end
            
            %Store F
            F_arr(num(k),k) = F;
        
            %increment F if solution is still stable
            if stab_flag == 1
                F = F + Fdir * dF;
                %increment counter
                num(k) = num(k) + 1;
            end
            %if unstable but at large dF, shrink dF and try again
            if stab_flag == 0 && dF == s*dF_0
                dF = dF_0;
                stab_flag = 1;
                F = F - Fdir * dF * (2*s + 2);
                if F * Fdir < F_0 * Fdir
                    F = F_0 + Fdir * dF;
                    num(k) = 2;
                else
                    num(k) = num(k) - 2;
                    [E_t1,~,~,~] = calc_eigenvalue(alpha,F,L,N,uout(:,num(k)-1));
                    [E_t2,~,~,~] = calc_eigenvalue(alpha,F + Fdir * dF,L,N,uout(:,num(k)-1));
                    while E_t1 > 0 || E_t2 > 0 || (E_t2 - E_t1) * eig_diff(num(k)-1,k) < 0
                        F = F - Fdir * dF * s;
                        num(k) = num(k) - 1;
                        [E_t1,~,~,~] = calc_eigenvalue(alpha,F,L,N,uout(:,num(k)-1));
                        [E_t2,~,~,~] = calc_eigenvalue(alpha,F + Fdir * dF,L,N,uout(:,num(k)-1));
                        if num(k) < 1
                            error("initial conditions too close to stability boundary");
                        end
                    end
                
                end
            end
            
        end
        
        %make sure enough stable point found to form F point set
        if num(k) < 3
            error("initial conditions too close to stability boundary");
        end

        %save last stable solution
        sol_init(:,k) = uout(:,num(k)-1);
        
        if k == 1
            uout_1 = uout;
        else
            uout_2 = uout;
        end
        
        %reset loop variables
        alpha = alpha + adir * da;
        F = F_0;
        uout = zeros(2*N,1);
        stab_flag = 1;
        dF = s*dF;
    
    end
    
    %move end number back until reached largest real eigenvalue part,
    num = num - 1;
    for k = 1:2
        while eig_diff(num(k),k) * eig_diff(num(k)-1,k) < 0 || E(num(k),k) > 0
            num(k) = num(k) - 1;
        end
        if k == 1
            sol_init(:,k) = uout_1(:,num(k));
        else
            sol_init(:,k) = uout_2(:,num(k));
        end
    end
    
    bif = bif(num(1),1);
    if bif == 1
        num = num + 1;
    end
    
    F1 = [F_arr(num(1),1)-2*Fdir*dF_0 F_arr(num(1),1)-1*Fdir*dF_0 F_arr(num(1),1)];
    F2 = [F_arr(num(2),2)-2*Fdir*dF_0 F_arr(num(2),2)-1*Fdir*dF_0 F_arr(num(2),2)];
    
    %create alpha and F vectors
    alpha_init = [alpha_0 alpha_0 + adir * da;];
    F_init = [round(F1,7); round(F2,7)];
    
    %call interpolation/extrapolate function on initial point set
    [F_bound_curr,err,F_init(1,:),F_sol(:,1),e_val(1,:),~] = boundary_locate(1,alpha_init(1),dF,F_init(1,:),Fdir,L,N,sol_init(:,1),bif,1);
    if F_bound_curr == -1
        error("unable to calculate initial boundary points. please choose different initial conditions")
    else
        fprintf('alpha = %15.10f calculated\n',alpha_init(1));
    end
    
    [F_bound_next,err,F_init(2,:),F_sol(:,2),e_val(2,:),~] = boundary_locate(1,alpha_init(2),dF,F_init(2,:),Fdir,L,N,sol_init(:,2),bif,1);
    if F_bound_next == -1
        error("unable to calculate initial boundary points. please choose different initial conditions")
    else
        fprintf('alpha = %15.10f calculated\n',alpha_init(2));
    end
    
    %rest of set up
    F_slope = (F_bound_next - F_bound_curr) / da;
    F_points = F_init(2,:) + F_slope * da;
    %array to store boundary points, starts equal to only calculated point
    F_bound = [F_bound_curr F_bound_next];
    alpha_bound = alpha_init;
    %solutions to return
    sol = F_sol;
    
end
