%% Model G_1 for diffusive dynamics
%% Semi-infinite BC
%% for source-detector pairs on the same plane (z=0) only!
%%
function g=G1Fitx(D,r,taustmp,muspo,muao,k0,ze)
	z0=1.0;
	r1=sqrt(r*r+z0^2./(muspo.*muspo));
	r2=sqrt(r*r+(z0+2*ze).^2./(muspo.*muspo));
   g=exp(-sqrt(3*muspo.*muao+6*k0.*k0.*muspo.*muspo.*D.*taustmp).*r1)./r1 ...
		-exp(-sqrt(3*muspo.*muao+6*k0.*k0.*muspo.*muspo.*D.*taustmp).*r2)./r2;
    g=g./(exp(-sqrt(3*muspo.*muao).*r1)./r1-exp(-sqrt(3*muspo.*muao).*r2)./r2);
%% revised: 2/25/10 EB