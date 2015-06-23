function X = muDPFfit(mu2ratio, mua0, musp0, rho, phi, theta)
    r = theta/phi;
    
    mua2 = mu2ratio(1)*mua0;
    musp2 = mu2ratio(2)*musp0;
    mua1 = (1-mu2ratio(1))*mua0;
    musp1 = (1-mu2ratio(2))*musp0;
    
    DPF0=3.*musp0./2./(rho*sqrt(3*mua0.*musp0)+1);    
    DPF1=3.*musp1./2./(rho*sqrt(3*mua1.*musp1)+1);
    DPF2=3.*musp2./2./(rho*sqrt(3*mua2.*musp2)+1);
    
    X = abs((DPF1/DPF0)*(1-r)*(mua1+musp1) + (DPF2/DPF0)*r*(mua2+musp2) - (mua0+musp0)*DPF0);    
end