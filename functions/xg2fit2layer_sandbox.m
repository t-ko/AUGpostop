%% Chi^2 fn for fitting G_2 (G2Fitx) (TK 10/29/14)
function X2=xg2fit2layer_sandbox(Db0,n0,R,mua,musp,lambda,ell,rho,taustmp,beta,g2,sigma)

%Forward solve for g1 with these values	
g1fit=g1fit2layer(Db0,n0,R,mua,musp,lambda,ell,rho,taustmp);
g1 = sqrt(abs((g2-1)./beta));

%Calculate difference between given solution and actual data
g1fit = reshape(g1fit,size(g1));
figure(9);semilogx(taustmp, g1fit); hold on; semilogx(taustmp, g1,'black');hold off
X2=norm((g1fit-g1).*sigma(1:length(g1)))