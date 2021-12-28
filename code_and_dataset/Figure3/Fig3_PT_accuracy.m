clear
close all
open Temperature_comparison2.fig

set(gca,'position',[0.1300 0.1100 0.35 0.8150])
colorbar
title('')
colorbar
text(-0.0099,-40,'B.','fontweight','b','fontsize',20)
cmocean('-grey')
grid on
grid minor
a1 = gca;

open Pressure_comparison3.fig

set(gca,'position',[0.1300 0.1100 0.35 0.8150])
colorbar
title('')
colorbar
text(-3.99,-40,'A.','fontweight','b','fontsize',20)
cmocean('-grey')
grid on
grid minor

a2 = copyobj(a1,gcf);
set(a2,'position',[0.5500 0.1100 0.4 0.8150])
axes(a2)
colorbar

set(gcf,'color','w')
export_fig -r200 PT_accuracy.png

fsize(16)
