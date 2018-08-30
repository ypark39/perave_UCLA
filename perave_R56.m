   for ip = 1:npart
        thetap(ip) = Ctot(ip+npart+(islice-1)*6*npart);
        phiabs_fin(ip) = phiabs(ip) + 2*pi*R56buncher/lambda*(gamma(ip)-gamma_avg)/gamma_avg+phaseshift;        
