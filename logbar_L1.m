function [y,container] = logbar_L1(y_init, phi, phi_t, z, nu, eps, tol_logbar,...
                               tol_newton, newton_iterlimit, tol_conjgrad, ...
                               conjgrad_iterlimit,newton_step, alpha, beta)
    tot_newton_iters = 0; % total newton iteration counter                        
    y=y_init; % initial point                       
    N = length(y_init);
    % disturbances
    no1 = 1.03;
    no2 = 0.015;
    % l1 norm
    w = no1*abs(y_init) + no2*max(abs(y_init));
    m = 2*length(y_init)+1;
    % tau chosen prudently, computed such that duality gap is roughly equal 
    % to <c^0, y> (original norm)after logbar iteration 1
    tau = (m)/sum(abs(y_init));
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
        w_newton = w;
        y_newton = y; % newton subroutine initial points
        res = phi(y_newton) - z; % residual
        fw_1 = y_newton - w_newton;
        fw_2 = -y_newton - w_newton;
        feps = (res'*res - eps^2)*0.5;
        f = sum(w_newton) - (1/tau)*(sum(log(-fw_1)) + sum(log(-fw_2)) + log(-feps));
        while (~completed)
            % follows equations in report
            sigma_11 = 1./fw_1.^2 + 1./fw_2.^2;
            sigma_12 = -1./fw_1.^2 + 1./fw_2.^2;
            sigma_Y = sigma_11 - sigma_12.^2./sigma_11;
            % bottom of Newton system
            nwgw = -tau - 1./fw_1 - 1./fw_2;
            phi_trans = phi_t(res);
            % top of Newton system
            nwgy = 1./fw_1 - 1./fw_2 + 1/feps*phi_trans;
            nablaf = -(1/tau)*[nwgy; nwgw];
            % xi 1 problem
            xi_1 = nwgy - sigma_12./sigma_11.*nwgw;
            % helper function to compute Hessian_11         
            Hessian_11_func = @(o) sigma_Y.*o - (1/feps)*phi_t(phi(o)) ...
                           + 1/feps^2*(phi_trans'*o)*phi_trans;
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
                q_cg = Hessian_11_func(d_cg); 
                alpha_cg = delta_cg/(d_cg'*q_cg);                
                z_cg = z_cg + alpha_cg*d_cg;                  
                if (mod(cg_iters+1,res_every) == 0)
                    res_cg = xi_1 - Hessian_11_func(z_cg);
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
            phi_dy = phi(dy);
            dw = (1./sigma_11).*nwgw - (sigma_12./sigma_11).*dy;              
            y_newton_next = y_newton + step*dy;
            w_newton_next = w_newton + step*dw;                        
            res_next = res + step*phi_dy;             
            while ( (res_next'*res_next > eps^2)|| ...
                            (max(abs(y_newton_next)-w_newton_next) > 0))
                iterations_cone = 1+iterations_cone;
                step = beta*step;                  
                res_next = res + step*phi_dy; 
                w_newton_next = w_newton + step*dw;
                y_newton_next = y_newton + step*dy; 
            end  
            % after centering, backtrack
            iters_backtracking = 0; % back tracking line search iterations
            fw_1_next = y_newton_next - w_newton_next;
            fw_2_next = -y_newton_next - w_newton_next;
            feps_next = (res_next'*res_next - eps^2)*0.5;            
            f_next = sum(w_newton_next) - (1/tau)*(sum(log(-fw_1_next)) + ...
                         sum(log(-fw_2_next)) + log(-feps_next));
            f_l = f + alpha*step*(nablaf'*[dy; dw]);                        
            while (f_l < f_next)
                iters_backtracking = iters_backtracking + 1;
                step = beta*step;
                res_next = res + step*phi_dy; 
                w_newton_next = w_newton + step*dw;
                y_newton_next = y_newton + step*dy; 
                f_l = f + alpha*step*(nablaf'*[dy; dw]); 
                fw_1_next = y_newton_next - w_newton_next;
                fw_2_next = -y_newton_next - w_newton_next;
                feps_next = (res_next'*res_next - eps^2)*0.5;
                f_next = sum(w_newton_next) - (1/tau)*(sum(log(-fw_1_next)) + ... 
                         sum(log(-fw_2_next)) + log(-feps_next));                                           
            end        
            newton_iters = 1+newton_iters;
            % updates for upcoming Newton iteration
            res = res_next;
            y_newton = y_newton_next; 
            w_newton = w_newton_next;                 
            f = f_next;
            fw_1 = fw_1_next; 
            fw_2 = fw_2_next;
            feps = feps_next; 
            lambda_2 = -(nablaf'*[dy; dw]);        
            % dummy variable 1 if Newton subroutine due to terminate
            % Newton subroutine terminates when Newton decrement less than or equal
            % to Newton tolerance or number of iters reaches
            % newton_iterlimit
            completed = (lambda_2/2 < tol_newton) | (newton_iters >= newton_iterlimit);    
            fprintf('LB iteration: %d. Newton iteration: %d. Duality gap: %f \n', count1, newton_iters, m/tau)          
        end
        tau = nu*tau;
        w = w_newton; 
        y = y_newton;
        tot_newton_iters = newton_iters + tot_newton_iters;
        container(count1, 1) = tot_newton_iters; % cumulative newton iterations
        container(count1, 2) =  m/tau; % duality gap
        fprintf('LB iteration %d undertaken. Total number of Newton iterations so far: %d \n', count1, tot_newton_iters)                        
    end                                            
end
                           