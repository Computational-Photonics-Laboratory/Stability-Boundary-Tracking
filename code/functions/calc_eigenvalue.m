function [E,uout,bif,res_flag] = calc_eigenvalue(alpha,F,L,N,sol)

    %First part of generating wave profile, takes place outside of loop
    dL = L/N;
    x = (-N/2:N/2-1)'*dL;
    % Fourier spectral differentiation d_xx
    %calculates second spatial derivative of del phi and del phi conj
    [~, D2_fourdif] = fourdif(N, 2);
    D2 = (2*pi/L)^2*D2_fourdif;

    % call fsolve
    options = optimset('Display','off','Jacobian','on','MaxIter',50);
    [uout,fval] = fsolve(@(u) LLPR_rhs(u,alpha,F,N,D2),sol,options);
    %calculate residual
    res = sum(sum(fval.*fval));
    
    %check residual - if above threshold, solution does not exist
    thresh = 1e-9; res_flag = 0;
    if res > thresh
        res_flag = 1;
    end
    
    %extract real and imaginary parts of wave then combine to get wave profile
    urout = uout(1:N);
    uiout = uout(N+1:N*2);
    u_ans = urout + 1j*uiout;
    u_abs = abs(u_ans).^2;
    
    %calculate range of intensity
    u_max = max(u_abs); u_min = min(u_abs);
    diff = u_max - u_min;
    
    %calculate eigenvalues
    phi_0 = u_ans .* eye(N);
    %create L matrix
    L11 = i*(D2+2*abs(phi_0)^2-alpha*eye(N))-eye(N);
    L12 = i*phi_0^2;
    L21 = conj(L12);
    L22 = conj(L11);
    L = [L11 L12; L21 L22];
    [~,D] = eig(L);
    D = diag(D);
    %take real e-values of L
    e_r = real(D);
    e_i = imag(D);

    %filter out values under threshold and determine largest real part
    thresh = 1e-7;
    num_max = 32;
    [er_max,ind] = maxk(e_r,num_max);
    max_arr = [0 0 0]; pos_arr = [0 0 0];
    m = 1;
    filter_count = 0;
    for n = 1:num_max
        if abs(er_max(n)) > thresh && max_arr(m) == 0
            max_arr(m) = er_max(n);
            pos_arr(m) = ind(n);
            m = m + 1;
        end
        %break if top three already filled
        if m > 3
            break;
        end
    end
    
    if n > 3
        filter_count = n - 3;
    end
    
    %check if the largest is the 0 value that has error above thresh
    %by seeing if the next largest is negative or not
    %for cases where the 0 error can exceed the threshold and hide 
    %the critical eigenvalue (pair)
    if filter_count == 0 && max_arr(1) > 0 && max_arr(2) < 0
        r_max = max_arr(2);
        pos_max = pos_arr(2);
    elseif filter_count == 0 && max_arr(1) < 0
        r_max = max_arr(2);
        pos_max = pos_arr(2);
    else
        r_max = max_arr(1);
        pos_max = pos_arr(1);
    end
    
    %set bifurcation type based on imaginary value of eigenvalue
    %with largest real part
    %if above thresh - complex conjugate pairs - hopf
    %below thresh - single values at Im(lambda) = 0 - saddle node
    if abs(e_i(pos_max)) > thresh
        bif = 1;
    else
        bif = 2;
    end
        
    %assign return value
    E = r_max;
    
    thresh = 1e-2;
    if bif == 2 && E > 0
        E = 1;
    %set E to 1 if residual is too large over a saddle node bifurcation
    elseif bif == 2 && res_flag == 1
        E = 2;
    %check if amplitude of intensity is above a threshold - collapse to cw
    elseif bif == 2 && diff < thresh
        E = 3;
    end
    
end
