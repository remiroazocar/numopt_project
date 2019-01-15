function z = phi_fourier_half_plane(x, om)
    N = round(sqrt(length(x)));
    z_c = 1/N*fft2(reshape(x,N,N));
    % returns real and im components of 2D fast fourier transform on upper half plane 
    % of transform domain.    
    z = [z_c(1,1); real(z_c(om))*sqrt(2); imag(z_c(om))*sqrt(2)];
end




