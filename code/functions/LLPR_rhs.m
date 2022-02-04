function [F,J] = LLPR_rhs(u,theta,h,N,D2)

ur = u(1:N);
ui = u(N+1:N*2);

F1 = D2*ur - ui - theta*ur + (ur.^2 + ui.^2).*ur;
F2 = D2*ui + ur - theta*ui + (ur.^2 + ui.^2).*ui - h;

F = [F1; F2];

if nargout > 1  
    e = ones(N,1);
    e = sparse(1:N, 1:N, e, N, N);      % identity matrix
    
    J = [D2 - theta*e + diag(3*ur.^2+ui.^2),...
           -e + diag(2*ur.*ui);...
           e + diag(2*ur.*ui),...
           D2 - theta*e + diag(ur.^2+3*ui.^2)];  
end

end
