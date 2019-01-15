function x = phi_t_fourier_half_plane(z, om, N)
    func = zeros(N,N); func(1,1)=z(1); k = length(z);
    func(om) = sqrt(2)*(z(2:(k+1)/2) + 1i*z((k+3)/2:k));
    x = reshape(real(N*ifft2(func)), N*N, 1);
    % returns real and im components of 2D fast fourier transform on bottom half plane of 
    % transform domain.
end    
    