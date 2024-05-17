function smeared_phase = smear_phase(phase,psf_um,gridspacing_um,smear_dim)

if psf_um > 0
    
    nr_gridpoints = size(phase,smear_dim);
    
    % smear function: gaussian
    smeargrid=(-(nr_gridpoints-1)/2:1:(nr_gridpoints-1)/2)*gridspacing_um;
    smearfunc=exp(-smeargrid.^2/(2*psf_um^2)); % gaussian at gridpoints, not normalized
    smearfunc=smearfunc/sum(smearfunc);
    
    smearfunc_fft = fft(smearfunc);
    smearfunc_fft_cloned = array_of_vectors(smearfunc_fft,size(phase),smear_dim);
    
    phase_fft = fft(phase,[],smear_dim);
    
    smeared_phase=ifft(smearfunc_fft_cloned.*phase_fft,[],smear_dim);
    smeared_phase=fftshift(smeared_phase,smear_dim);
    
else
    
    smeared_phase = phase;
    
end

end