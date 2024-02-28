%script explain distance
ccc
N=1:200;
M=N;
q=.0005;
r=.5;

%  pnode=0;
%  player=0;


 % %uniform multiplex
 pnode_mp = (N.*(N-1)-2)./(M.*N.*(N-1)-2); 
 player_mp = 2*(M-1)./(M.*N.*(N-1)-2);  
 
% %uniform multilayer
 pnode_ml = (N.*N-2)./((M.*N).^2-2) ; 
 player_ml= 2*(M.^2-1)./((M.*N).^2-2) ;

%  pnode_ml=0;
%  pnode_mp=0;
%  player_ml=0;
%  player_mp=0;

ptel_mp = 1-pnode_mp-player_mp;
ptel_ml = 1-pnode_ml-player_ml;

% %multilayer
 dx_ml=r.*(pnode_ml+ptel_ml).*N.*M'.*sqrt(q.*(1-q));
 dy_ml=r.*(player_ml+ptel_ml).*N.*M'.*sqrt(q.*(1-q));

% %multiplex layer ??? still not clear
 dx_mp=r.*(pnode_mp+ptel_mp).*M.*sqrt(N)'.*sqrt(q.*(1-q));
 dy_mp=r.*(player_mp+ptel_mp).*M.*sqrt(N)'.*sqrt(2*q.*(1-q));

figure
subplot(2,2,1); surf(dx_ml); colormap(jet); colorbar ; axis square tight
caxis([0 max(max(max(dx_ml)),max(max(dy_ml)))]);
shading interp; view(2);
xlabel('N');
ylabel('M');
title('d_X ML');
subplot(2,2,2); surf(dy_ml); colormap(jet); colorbar ; axis square tight
caxis([0 max(max(max(dx_ml)),max(max(dy_ml)))]);
shading interp; view(2);
xlabel('N');
ylabel('M');
title('d_Y ML');
subplot(2,2,3); surf(dx_mp); colormap(jet); colorbar ; axis square tight
caxis([0 max(max(max(dx_mp)),max(max(dy_mp)))]);
shading interp; view(2);
xlabel('N');
ylabel('M');
title('d_X MP');
subplot(2,2,4); surf(dy_mp); colormap(jet); colorbar ; axis square tight
caxis([0 max(max(max(dx_mp)),max(max(dy_mp)))]);
shading interp; view(2); 
xlabel('N');
ylabel('M');
title('d_Y MP');
%%%
dx_ml_f=diag(fliplr(dx_ml));
dy_ml_f=diag(fliplr(dy_ml));
dx_mp_f=diag(fliplr(dx_mp));
dy_mp_f=diag(fliplr(dy_mp));

figure 
yyaxis right
plot([max(N)/3 max(N)/3],[0 max(dy_mp_f)],'--k')
hold on
plot(N-fliplr(M),fliplr(dx_mp_f'),'-','linewidth',1.5,'color',[189/255 215/255 238/255]);
hold on
yyaxis right
plot(N-fliplr(M),fliplr(dy_mp_f'),'-','linewidth',1.5,'color',[248/255 203/255, 173/255]);
ylim([0 max(dx_ml_f)/2])
ylabel('Multiplex distance')

yyaxis left
plot(N-fliplr(M),fliplr(dx_ml_f'),'-','linewidth',1.5,'color',[.75 .075 .75]);
hold on

yyaxis left
plot([0 0],[0 max(dy_ml_f)],'--k')
hold on
plot(N-fliplr(M),fliplr(dy_ml_f'),'-k','linewidth',1.5,'color',[.25 .25 .25]);
ylim([0 max(dx_ml_f)])
xlabel('Nodes-Layers (N-M)');
ylabel('Multilayer distance')
%legend('d_X','d_Y');

ax=gca;
yyaxis left
ax.YColor=[.15 .15 .15];

yyaxis right
ax.YColor=[.75 .15 .75];

figure(3)
loglog(dx_ml(end,:),'-b');
hold on
loglog(dy_ml(end,:),'-r');
xlabel('M')
ylabel('Distances')
legend('d_X','d_Y')
% figure(4)
% 
% aux=(1-pnode_ml)./(1-player_ml);
% plot(aux)

%axis square

% figure
% plot(n,dx');
% hold on
% dy2=fliplr(dy);
% plot(n,dy2','r');
%  hold on
% plot(n,dy2'-dx');
%figure 
%subplot(1,2,1); imagesc(dx); colormap(jet); colorbar ; axis square
%subplot(1,2,2); imagesc(dy); colormap(jet); colorbar ; axis square

%figure
 %plot(dx); colormap(jet); colorbar ; axis square
%imagesc(dy-dx); colormap(jet); colorbar ; axis square
