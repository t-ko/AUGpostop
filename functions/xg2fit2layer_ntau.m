%% Chi^2 fn for fitting G_2 (G2Fitx) (TK 10/29/14)
function X2=xg2fit2layer_ntau(Db0,n0,R,mua,musp,lambda,ell,rho,taus,beta,g2tmp,startcorr,tmpf,n_tau)
incr = floor(length(taus(startcorr:tmpf))./(n_tau-1));

%Forward solve for g1 with these values	
g1fit=g1fit2layer(Db0,n0,R,mua,musp,lambda,ell,rho,taus([1 startcorr:incr:tmpf end]));
g1fit = g1fit(2:end-1);
g1 = sqrt(abs((g2tmp(1:incr:end)-1)./beta));

%Calculate difference between given solution and actual data
g1fit = reshape(g1fit,size(g1));
%figure(10);semilogx(taus(startcorr:incr:tmpf), g1fit); hold on; semilogx(taus(startcorr:incr:tmpf), g1,'black');hold off
% X2=norm((g1fit-g1)./sigmatmp(1:incr:end));
X2=norm((g1fit-g1));
