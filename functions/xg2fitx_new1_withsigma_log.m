%% Chi^2 fn for fitting G_2 (G2Fitx)
%% t and G2 are column vectors
%%
%% x0(1)=D_B , x0(2)=beta
%%revised by EB 6/5/09
%%log-weighting added by TK 08/20/14
function g=xg2fitx_new1_withsigma_log(Db,beta,r,taustmp,muspo,muao,k0,ze,G2,N,sigma)
    x0=[Db beta];
	f=g2fitx(x0,r,taustmp,muspo,muao,k0,ze);
    
	%g=norm(f-G2(1:N));
    %g=norm((f-G2(1:N)')./sigma(1:N)).*log(taustmp);
    g=norm((f-G2(1:N)')./sigma(1:N).*log(taustmp));