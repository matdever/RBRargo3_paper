clear

% load all the data using the WHOI routine to import WHOI-decoded data
%data = WHOI_float_import('/Users/mdever/Dropbox (RBR)/1RBR Science/Argo/data/WMO4903275/');
% ~10 cm/s --> profile 1-36 (32 is nice)
% ~20 cm/s --> profile 37-74 (71 is nice)
% ~5 cm/s --> profile 81-161 (81 is nice)
% ~15 cm/s --> profile 162-191 (81 is nice)
% ~3 cm/s --> profile 192-212 (81 is nice)
% ~15 cm/s --> profile 213-223 (81 is nice)
% ~10 cm/s --> profile 225-333 (81 is nice)


%%
clear
load('4903275.mat')
counter = 1;
for ii = [ 135 306 189 44]%206
    ii
    
    % Extract variables
    try
        Tmeas = data{ii}.prof.temp;
    catch
        continue
    end
    Pmeas = data{ii}.prof.pres;
    Smeas = data{ii}.prof.psal;
    Tcond = data{ii}.prof.cttemp;
    
    ind = find(data{ii}.rise.pres<200,1,'first');
    ind2 = find(data{ii}.rise.pres<5,1,'first');
    etime = interp1(data{ii}.rise.pres(ind:ind2),data{ii}.rise.sec(ind:ind2),Pmeas);
    dpdt = cat(2,NaN,(Pmeas(3:end) - Pmeas(1:end-2))./(etime(3:end)-etime(1:end-2)),NaN);
    
    % compute derived variables
    Cmeas = gsw_C_from_SP(Smeas,Tmeas,Pmeas);
    SAmeas = gsw_SA_from_SP(Smeas,Pmeas,data{ii}.gps(end).lon,data{ii}.gps(end).lat);
    CTmeas = gsw_CT_from_t(SAmeas,Tmeas,Pmeas);
    sigmeas = gsw_sigma0(SAmeas,CTmeas);
    %% apply dynamic corrections
    CTlag = -0.35;
    
    Tcor = NaN*Tmeas;
    good = find(~isnan(etime) & ~isnan(Tmeas));
    [~,IA,~] = unique(etime(good));
    Tcor(good(IA)) = interp1(etime(good(IA)), Tmeas(good(IA)), etime(good(IA))-CTlag);
    
    ref_vel = 0.101;
    ref_coef = 9.3e-3;
    ctcoef = 0.08*(-dpdt*100).^-0.89;
    Tlong = ctcoef.*(Tcond-Tmeas);
    
    ref_vel = 0.12;
    alpha = 0.53*(-dpdt*100).^-1.12;
    tau = 14.35*(-dpdt*100).^-0.24;
    fs = 1;
    fn = fs/2;
    
    a = 4*fn.*alpha.*tau./(1+4*fn.*tau);
    b = 1-2.*a./alpha;
    
    Tshort = NaN*Tmeas;
    % prep Tcor --> interpolate onto a 1hz grid and over NaNs
    good = find(~isnan(Tcor) & ~isnan(etime) & ~isnan(a) & ~isnan(b));
    [~,IA,~] = unique(etime(good));
    etimesyn = min(etime(good(IA))):fs:max(etime(good(IA)));
    asyn = interp1(etime(good(IA)),a(good(IA)),etimesyn);
    bsyn = interp1(etime(good(IA)),b(good(IA)),etimesyn);
    Tsyn = interp1(etime(good(IA)),Tcor(good(IA)),etimesyn);
    Tshortsyn = NaN*Tsyn;
    
    for tt = 1:length(Tshortsyn)
        if tt == 1
            Tshortsyn(tt) = 0;
        else
            Tshortsyn(tt) = -bsyn(tt).*Tshortsyn(tt-1) + asyn(tt).*(Tsyn(tt)-Tsyn(tt-1));
        end
    end
    % re-interpolates on orgiinal time axis
    Tshort(good(IA)) = interp1(etimesyn,Tshortsyn,etime(good(IA)));
    %clear *syn
    
    Scor = gsw_SP_from_C(Cmeas,Tcor+Tlong-Tshort,Pmeas);
    SAcor = gsw_SA_from_SP(Scor,Pmeas,data{ii}.gps(end).lon,data{ii}.gps(end).lat(end));
    CTcor = gsw_CT_from_t(SAcor,Tcor,Pmeas);
    sigcor = gsw_sigma0(SAcor,CTcor);
    
    %% FIGURE
    if counter == 1
        fullfigure
    end
    subplot(1,4,counter)
    plot(Smeas,Pmeas,'k','linewidth',2)
    grid on
    axis ij
    fsize(28)
    hold on
    plot(Scor,Pmeas,'r','linewidth',2)
    set(gca,'xaxislocation','top')
    xlabel('S_P [ ]')
    %xlim([34.9 35.01])
    if ii == 206
        legend('raw','corrected')
        ylabel('Pressure [dbar]')
        ylim([0 50])
        xlim([34.9 35.151])
        text(34.915,48.3,'3 cm/s','fontsize',25,'backgroundcolor','w','edgecolor','k')
        
    elseif ii == 135
        legend('raw','corrected')
        ylabel('Pressure [dbar]')        
        xlim([34.9 35.151])
        ylim([0 60])
        text(34.915,58,'5 cm/s','fontsize',25,'backgroundcolor','w','edgecolor','k')

    elseif ii == 306
        ylim([0 60])
        xlim([34.9 35.151])
        text(34.915,58,'10 cm/s','fontsize',25,'backgroundcolor','w','edgecolor','k')
    
    elseif ii == 189
        xlim([34.9 35.151])
        ylim([0 45])
        text(34.915,43.5,'15 cm/s','fontsize',25,'backgroundcolor','w','edgecolor','k')

    elseif ii == 44
        xlim([34.9 35.151])
        ylim([0 45])
        text(34.915,43.5,'20 cm/s','fontsize',25,'backgroundcolor','w','edgecolor','k')
        
    end
    set(gca,'xtick',34.9:0.1:35.15)
    
    grid minor
    counter = counter + 1;
    set(gcf,'color','w')
    
end
export_fig('Fig13_Arcticfloat.png','-r200')
