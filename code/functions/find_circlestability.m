function [stab_indices,bif,sol_init,err] = find_circlestability(alpha,F,L,N,sol)

    sol_init = sol;
    err = 0;
    
    num_points = length(alpha);

    %First part of generating wave profile, takes place outside of loop
    dL = L/N;
    x = (-N/2:N/2-1)'*dL;
    % Fourier spectral differentiation d_xx
    %calculates second spatial derivative of del phi and del phi conj
    [~, D2_fourdif] = fourdif(N, 2);
    D2 = (2*pi/L)^2*D2_fourdif;
    options = optimset('Display','off','Jacobian','on','MaxIter',50);
    
    %stability flag - exit for loop upon finding unstable point,
    %assuming all remaining points are also unstable
    stab_flag = 1;
    
    %bifurcation type vector
    bif = ones(num_points,1);
    
    %loop counter
    p = 1;
    
    while p <= num_points && stab_flag == 1
        
        %get largest real eigenvalue and bifurcation type of current point
        if p == 1
            [E,uout(:,p),bif(p),~] = calc_eigenvalue(alpha(p),F(p),L,N,sol);
        else
            [E,uout(:,p),bif(p),~] = calc_eigenvalue(alpha(p),F(p),L,N,uout(:,p-1));
        end
        
        %if largest real eigenvalue part above 0, set stability flag to 0
        if E > 0
            stab_flag = 0;
        end
        
        %set respective array index to stable if point is stable
        if stab_flag == 0
            if bif(p) == 2
                %saddle node of transcritical bifurcation
                if p >= 4
                    stab_indices = [p-3 p-2 p-1];
                    bif = 2;
                    if p > 4
                        sol_init = uout(:,p-4);
                    end
                else
                    stab_indices = [0 0 0];
                    bif = 0;
                    err = 1;
                end
            else
                %hopf bifurcation
                if p >= 3
                    stab_indices = [p-2 p-1 p];
                    bif = 1;
                    if p > 3
                        sol_init = uout(:,p-3);
                    end
                else
                    stab_indices = [0 0 0];
                    bif = 0;
                    err = 1;
                end
            end
        end
        
        %increment counter
        p = p + 1;
       
    end

    %throw error if no unstable points detected
    if stab_flag == 1
        bif = 0;
        err = 1;
        stab_indices = [0 0 0];
    end

end