function [y, container]= logbar_TV(y_init, phi, phi_t, z, nu, eps, tol_logbar,...
                                   tol_newton, newton_iterlimit, tol_conjgrad, ...
                                   conjgrad_iterlimit,newton_step, alpha, beta)

    tot_newton_iters = 0; % total newton iteration counter                        
    y=y_init; % initial point 
    % compose D_{h;ij} and D_{v;ij} sparse matrices for Total-Variation                        
    N = length(y_init);    
    Dh = spdiags([reshape([-ones(round(sqrt(N)),round(sqrt(N))-1)...
                    zeros(round(sqrt(N)),1)],N,1) ...
                    reshape([zeros(round(sqrt(N)),1) ones(round(sqrt(N)),...
                    round(sqrt(N))-1)],N,1)], [0 round(sqrt(N))], N, N);
    Dv = spdiags([reshape([-ones(round(sqrt(N))-1,round(sqrt(N))); 
                    zeros(1,round(sqrt(N)))],N,1), ...
                    reshape([zeros(1,round(sqrt(N))); ones(round(sqrt(N))-1, ...
                    round(sqrt(N)))],N,1)], [0 1], N, N);
    Dhy = Dh*y;
    Dvy = Dv*y;
    % disturbances
    no1 = 1.03;
    no2 = 0.015;
    Droot = sqrt(Dvy.^2 + Dhy.^2);
    % TV seminorm
    u = no1*Droot + no2*max(Droot);
    m = length(y_init)+1;
    % tau chosen prudently, computed such that duality gap is roughly equal 
    % to <c^0, y> after logbar iteration 1
    tau = m/sum(Droot);
    % total log barrier iterations required calculated in advance. 
    tot_logbar_iters = (log(m)-log(tau)--log(tol_logbar))/log(nu);
    % result container with cumulative newton iterations and duality gap
    % for each log barrier iteration. 
    container = zeros(ceil(tot_logbar_iters),2);
    fprintf('%d LB iterations  are required. \n', ceil(tot_logbar_iters))
    for count1 = 1:ceil(tot_logbar_iters)                        
        % newton subroutine commences (centering: newton method with back tracking
        % line search along search dir        
        completed = 0;
        newton_iters = 0;                        
        u_newton = u;
        y_newton = y; % newton subroutine initial points
        res = phi(y_newton) - z; % residual
        % reset D_{h;ij} and D_{v;ij} sparse matrices for Total-Variation                  
        Dh = spdiags([reshape([-ones(round(sqrt(N)),round(sqrt(N))-1)...
                        zeros(round(sqrt(N)),1)],N,1) ...
                        reshape([zeros(round(sqrt(N)),1) ones(round(sqrt(N)),...
                        round(sqrt(N))-1)],N,1)], [0 round(sqrt(N))], N, N);
        Dv = spdiags([reshape([-ones(round(sqrt(N))-1,round(sqrt(N))); 
                        zeros(1,round(sqrt(N)))],N,1), ...
                        reshape([zeros(1,round(sqrt(N))); ones(round(sqrt(N))-1, ...
                        round(sqrt(N)))],N,1)], [0 1], N, N);
        % follows equations in report
        Dvy = Dv*y_newton;                 
        Dhy = Dh*y_newton;          
        feps = (res'*res - eps^2)*0.5;
        fu = (Dhy.^2 - u_newton.^2 + (Dvy).^2)*0.5;        
        f = sum(u_newton) - (1/tau)*(sum(log(-fu)) + log(-feps));
        while (~completed)
            sigma_12 = -u_newton./fu.^2;
            sigma_22 = 1./fu + (u_newton.^2)./(fu.^2);            
            sigma_Z = 1./fu.^2 - (sigma_12.^2)./sigma_22;                                    
            % bottom of newton system
            nugu = -tau - u_newton./fu;
            phi_trans = phi_t(res); 
            % top of newton system
            nugy = Dh'*((1./fu).*Dhy) + Dv'*((1./fu).*Dvy) + 1/feps*phi_trans;            
            nablaf = -(1/tau)*[nugy; nugu];            
            % xi 1 problem
            xi_1 = nugy - Dh'*(Dhy.*(sigma_12./sigma_22).*nugu) - Dv'*(Dvy.*(sigma_12./sigma_22).*nugu);
            % helper function to compute Hessian_11         
            Hessian_11_fun = @(o) Hessian_11(o, phi, phi_t, Dh, Dv, Dhy, Dvy, sigma_Z, fu, feps, phi_trans);
            % start of conjugate gradient subroutine
            % conjugate gradient to solve symmetric positive definite
            % system Hz = g where g is xi_1, z is z_cg, H is Hessian_11
            res_every = 50;
            cg_iters = 0; 
            res_cg = xi_1;
            d_cg = res_cg;
            delta_0_cg = xi_1'*xi_1;
            delta_cg = res_cg'*res_cg;
            z_cg = zeros(length(xi_1),1);
            top_z_cg = z_cg;
            top_res_cg = sqrt(delta_cg/delta_0_cg);
            % cg subroutine terminates when ||Hz - g||_2/||g||_2 less than 
            % tol_conjgrad or conjgrad_iterlimit iterations reached
            while ((delta_cg> tol_conjgrad^2*delta_0_cg) && ...
                                (cg_iters < conjgrad_iterlimit))
                q_cg = Hessian_11_fun(d_cg); 
                alpha_cg = delta_cg/(d_cg'*q_cg);                
                z_cg = z_cg + alpha_cg*d_cg;                  
                if (mod(cg_iters+1,res_every) == 0)
                    res_cg = xi_1 - Hessian_11_fun(z_cg);
                else
                    res_cg = res_cg - alpha_cg*q_cg;
                end
                cg_iters = cg_iters+ 1;
                prev_delta_cg = delta_cg; 
                delta_cg = res_cg'*res_cg;
                beta_cg = delta_cg/prev_delta_cg;
                d_cg = res_cg + beta_cg*d_cg;
                if (sqrt(delta_cg/delta_0_cg) < top_res_cg)
                    top_res_cg = sqrt(delta_cg/delta_0_cg);
                    top_z_cg = z_cg;                    
                end    
                
            end
            z_cg = top_z_cg;
            % end of conjugate gradients subroutine
            iterations_cone = 0;
            step = newton_step; % initial newton step is 0.9
            dy = z_cg;
            Dvdy = Dv*dy;
            Dhdy = Dh*dy;  
            phi_dy = phi(dy);            
            du = (1./sigma_22).*(nugu - sigma_12.*(Dhy.*Dhdy + Dvy.*Dvdy));            
            y_newton_next = y_newton + step*dy;
            u_newton_next = u_newton + step*du;                        
            res_next = res + step*phi_dy;  
            Dhy_next = Dhy + step*Dhdy;  
            Dvy_next = Dvy + step*Dvdy;            
            while ((res_next'*res_next > eps^2) || ...
                        (max(sqrt(Dhy_next.^2+Dvy_next.^2) - u_newton_next) > 0))
                iterations_cone = 1+iterations_cone;
                step = beta*step;                  
                Dvy_next = Dvy + step*Dvdy;
                Dhy_next = Dhy + step*Dhdy;
                res_next = res + step*phi_dy; 
                u_newton_next = u_newton + step*du;
                y_newton_next = y_newton + step*dy; 
            end            
            % after centering, backtrack
            iters_backtracking = 0; % back tracking line search iterations
            fu_next = (Dhy_next.^2 + Dvy_next.^2 - u_newton_next.^2)*0.5;
            feps_next = (res_next'*res_next - eps^2)*0.5;            
            f_next = sum(u_newton_next) - (1/tau)*(sum(log(-fu_next)) + log(-feps_next));
            f_l = f + alpha*step*(nablaf'*[dy; du]);            
            while (f_l < f_next)
                iters_backtracking = iters_backtracking + 1;
                step = beta*step;
                res_next = res + step*phi_dy;  
                Dhy_next = Dhy + step*phi_dy;
                Dvy_next = Dvy + step*phi_dy;                
                u_newton_next = u_newton + step*du;
                y_newton_next = y_newton + step*dy;                  
                f_l = f + alpha*step*(nablaf'*[dy; du]); 
                fu_next = (Dhy_next.^2 + Dvy_next.^2 - u_newton_next.^2)*0.5;
                feps_next = (res_next'*res_next - eps^2)*0.5;
                f_next = sum(u_newton_next) - (1/tau)*(sum(log(-fu_next)) + log(-feps_next));                             
            end
            newton_iters = 1+newton_iters;
            % updates for upcoming Newton iteration
            Dvy = Dvy_next;
            Dhy = Dhy_next; 
            res = res_next;
            y_newton = y_newton_next; 
            u_newton = u_newton_next;                 
            f = f_next;
            fu = fu_next; 
            feps = feps_next; 
            lambda_2 = -(nablaf'*[dy; du]);
            % dummy variable 1 if Newton subroutine due to terminate
            % Newton subroutine terminates when Newton decrement less than or equal
            % to Newton tolerance or number of iters reaches
            % newton_iterlimit
            completed = (lambda_2/2 < tol_newton) | (newton_iters >= newton_iterlimit);    
            fprintf('LB iteration: %d. Newton iteration: %d. Duality gap: %f \n', count1, newton_iters, m/tau)
        end      
        tau = nu*tau;
        u = u_newton; 
        y = y_newton;
        tot_newton_iters = newton_iters + tot_newton_iters;
        container(count1, 1) = tot_newton_iters; % cumulative newton iterations
        container(count1, 2) =  m/tau; % duality gap
        fprintf('LB iteration %d undertaken. Total number of Newton iterations so far: %d \n', count1, tot_newton_iters)
    end                                                                                                                                       
end


function z=Hessian_11(y, phi, phi_t, Dh, Dv, Dhx, Dvx,sigma_Z, fu,feps,phi_trans)
    Dvy = Dv*y;
    Dhy = Dh*y;
    z = Dh'*((-1./fu + sigma_Z.*Dhx.^2).*Dhy + sigma_Z.*Dhx.*Dvx.*Dvy) + ...
        Dv'*((-1./fu + sigma_Z.*Dvx.^2).*Dvy + sigma_Z.*Dhx.*Dvx.*Dhy) - ...
        1/feps*phi_t(phi(y)) + 1/feps^2*(phi_trans'*y)*phi_trans;  

end



