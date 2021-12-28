clear

load Stability_analysis.mat
profile_number = repmat(profile_number,19,1);
%%
fullfigure
%subplot(2,1,1)
axes('position', [0.1300    0.5838    0.7750    0.30])
hold on
grid
scatter(profile_number(:), diffsal0(:), 35,'ok', 'filled')
XLIM = get(gca,'xlim');
patch([XLIM(1) XLIM(1) XLIM(2) XLIM(2) XLIM(1)],[-0.01 0.01 0.01 -0.01 -0.01],[77,196,77]/256,'facealpha',0.3)
scatter(profile_number(:), diffsal0(:), 35,'ok', 'filled')
ylabel('∆S - ∆S (t=0) [PSU]','position',[-7 -0.05 -1])
%xlabel('Profile number')
fsize(30)
h1 = scatter(profile_number(3,:), diffsal0(3,:), 35,'or', 'filled');
ylim([-0.05 0.05])
box on
legend(h1,'WMO690299')
set(gca,'xticklabel',{})

axes('position', [0.1300    0.425    0.7750    0.15])
hold on
grid
xlim([0 140])
%ylabel('∆S - ∆S (t=0) [PSU]')
xlabel('Profile number')
fsize(30)
h1 = scatter(profile_number(3,:), diffsal0(3,:), 35,'or', 'filled');
ylim([-0.40 -0.05])
box on

set(gcf,'color','w')
export_fig -r200 Fig6_stability.png
% %%
% blop = diffsal0;
% blop(3,:) = [];
% boxplot(blop)
% hold on
% patch([XLIM(1) XLIM(1) XLIM(2) XLIM(2) XLIM(1)],[-0.01 0.01 0.01 -0.01 -0.01],'g','facealpha',0.2)
