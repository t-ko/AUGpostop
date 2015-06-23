%% Chi^2 fn for fitting G_2 (G2Fitx)
%% t and G2 are column vectors
%%
%% x0(1)=D_B , x0(2)=beta
%%revised by EB 6/5/09
function g=xg2fitx_betaandDB_new1_withsigma(x0,r,taustmp,muspo,muao,k0,ze,G2,N,sigma)
    
	f=g2fitx(x0,r,taustmp,muspo,muao,k0,ze);
    G2 = G2(1:N);
    G2 = reshape(G2,size(f));
	%g=norm(f-G2);
    g=norm((f-G2(1:N))./sigma(1:N));
    
