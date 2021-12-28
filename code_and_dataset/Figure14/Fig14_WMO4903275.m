clear

% matrixi = NaN*ones(3000,410);
% matrixP = NaN*ones(3000,410);
% matrixDSlong = NaN*ones(3000,410);
% matrixDSshort = NaN*ones(3000,410);
% matrixDSall = NaN*ones(3000,410);
% matrixDSsbe = NaN*ones(3000,410);
%

counter = 1;

% Load the data
WMO = 4903275;
for ii = [27:63 65:83 85:333]%[81:161 27:36 225:333 162:191 213:223 37:63 65:74]%
    filename = ['datasets/WMO',num2str(WMO),'/11149_',num2str(ii,'%03.0f'),'.lis'];
    theii(counter) = ii;
    try
        [TS, binneddata] = RBR_lis2mat(filename);
        baddata = find(isnan(binneddata.P) | isnan(binneddata.T) | isnan(binneddata.S) |isnan(binneddata.Tcond) | binneddata.Tcond==-999);
        binneddata(baddata,:) = [];
        % create timestamps
        binneddata.etime = interp1(TS.P,TS.etime,binneddata.P);
        
        % sort data
        [binneddata.etime,I] = sort(binneddata.etime);
        binneddata.T = binneddata.T(I);
        binneddata.P = binneddata.P(I);
        binneddata.Tcond = binneddata.Tcond(I);
        binneddata.S = binneddata.S(I);
        
        % computed variables
        binneddata.C = gsw_C_from_SP(binneddata.S,binneddata.T,binneddata.P);

        % average profiling speed
        dpdt(counter) = -nanmean(diff(binneddata.P(binneddata.P<100))./diff(binneddata.etime(binneddata.P<100))*100);
        
        ctcoef = 0.07*dpdt(counter).^(-0.86) + 7.5e-9;
        
        Tlong = ctcoef*(binneddata.Tcond-binneddata.T);
        binneddata.Scorlong = gsw_SP_from_C(binneddata.C,binneddata.T+Tlong,binneddata.P);

        
%         alpha = 0.041;
%         tau = 8.11;
         alpha = 0.5*dpdt(counter).^(-1.09) + 1.6e-7;
         tau = 13.46*dpdt(counter).^(-0.22) + 2.2e-5;
        % Interpolate needed variables onto a 1Hz time series
        gooddata = find(~isnan(binneddata.T) & ~isnan(binneddata.etime));
        gridded.etime = min(binneddata.etime):1:max(binneddata.etime)+1;
        [C,IA,~] = unique(binneddata.etime(gooddata));
        gridded.TEMP = interp1(binneddata.etime(gooddata(IA)),binneddata.T(gooddata(IA)),gridded.etime);
        
        % Samplig frequency
        fs = 1;
        % Nyquist frequency
        fn = fs/2;
        
        % Lueck and Picklo (1991)'s coefficients (a,b)
        a = 4*fn.*alpha.*tau./(1+4*fn.*tau);
        b = 1-2.*a./alpha;
        
        % Apply Lueck and Picklo (1991) recursive filter
        gridded.Tshort = zeros(size(gridded.TEMP));
        
        % Set first value to 0°C.
        gridded.Tshort(1) = 0;
        
        for tt = 2:length(gridded.TEMP)
            gridded.Tshort(tt) = -b*gridded.Tshort(tt-1) + a*(gridded.TEMP(tt)-gridded.TEMP(tt-1));
        end; clear tt
        
        
        % Go back and set the first temperature adjustment to the second, since the
        % float is probably rising as it takes its first samples
        gridded.Tshort(1)=gridded.Tshort(2);
        
        % Re-grid onto original time grid
        Tshort = NaN*binneddata.T;
        Tshort(gooddata) = interp1(gridded.etime,gridded.Tshort,binneddata.etime(gooddata));
        clear gridded a b fs fn gooddata
        
        binneddata.Scorshort = gsw_SP_from_C(binneddata.C,binneddata.T-Tshort,binneddata.P);
        binneddata.Scor = gsw_SP_from_C(binneddata.C,binneddata.T-Tshort+Tlong,binneddata.P);
        
        
        alpha = 0.141;
        tau = 6.68;
        % Interpolate needed variables onto a 1Hz time series
        gooddata = find(~isnan(binneddata.T) & ~isnan(binneddata.etime));
        gridded.etime = min(binneddata.etime):1:max(binneddata.etime)+1;
        [C,IA,~] = unique(binneddata.etime(gooddata));
        gridded.TEMP = interp1(binneddata.etime(gooddata(IA)),binneddata.T(gooddata(IA)),gridded.etime);
        
        % Sampling frequency
        fs = 1;
        % Nyquist frequency
        fn = fs/2;
        
        % Lueck and Picklo (1991)'s coefficients (a,b)
        a = 4*fn.*alpha.*tau./(1+4*fn.*tau);
        b = 1-2.*a./alpha;
        
        % Apply Lueck and Picklo (1991) recursive filter
        gridded.Tshort = zeros(size(gridded.TEMP));
        
        % Set first value to 0°C.
        gridded.Tshort(1) = 0;
        
        for tt = 2:length(gridded.TEMP)
            gridded.Tshort(tt) = -b*gridded.Tshort(tt-1) + a*(gridded.TEMP(tt)-gridded.TEMP(tt-1));
        end; clear tt
        
        % Go back and set the first temperature adjustment to the second, since the
        % float is probably rising as it takes its first samples
        gridded.Tshort(1)=gridded.Tshort(2);
        
        % Re-grid onto original time grid
        Tshort = NaN*binneddata.T;
        Tshort(gooddata) = interp1(gridded.etime,gridded.Tshort,binneddata.etime(gooddata));
        clear gridded a b fs fn gooddata
        
        binneddata.Scorsbe = gsw_SP_from_C(binneddata.C,binneddata.T-Tshort,binneddata.P);
        
        matrixi(:,counter) = counter;
        matrixprof(:,counter) = ii;
        matrixP(1:length(binneddata.P),counter) = binneddata.P;
        matrixDSlong(1:length(binneddata.P),counter) = (binneddata.S-binneddata.Scorlong);
        matrixDSshort(1:length(binneddata.P),counter) = (binneddata.S-binneddata.Scorshort);
        matrixDSall(1:length(binneddata.P),counter) = (binneddata.S-binneddata.Scor);
        matrixDSsbe(1:length(binneddata.P),counter) = (binneddata.S-binneddata.Scorsbe);
        
        binneddata2.P = 0:.05:200;
        [C,IA,~] = unique(binneddata.P);
        binneddata2.S = interp1(binneddata.P(IA),binneddata.S(IA),binneddata2.P);
        binneddata2.Scor = interp1(binneddata.P(IA),binneddata.Scor(IA),binneddata2.P);
        binneddata2.Scorlong = interp1(binneddata.P(IA),binneddata.Scorlong(IA),binneddata2.P);
        binneddata2.Scorshort = interp1(binneddata.P(IA),binneddata.Scorshort(IA),binneddata2.P);
        binneddata2.Scorsbe = interp1(binneddata.P(IA),binneddata.Scorsbe(IA),binneddata2.P);
        
        matrixi2(:,counter) = counter;
        matrixprof2(:,counter) = ii;
        matrixP2(:,counter) = binneddata2.P;
        matrixDSlong2(:,counter) = (binneddata2.S-binneddata2.Scorlong);
        matrixDSshort2(:,counter) = (binneddata2.S-binneddata2.Scorshort);
        matrixDSall2(:,counter) = (binneddata2.S-binneddata2.Scor);
        matrixDSsbe2(:,counter) = (binneddata2.S-binneddata2.Scorsbe);
        
        counter = counter +1;
        
    catch
        continue
    end
    
end
matrixi = repmat(matrixi,2247,1);
matrixi2 = repmat(matrixi(1,:),4001,1);

%%
% Sorat matrices by fall-rates
[dpdt,I] = sort(dpdt);
matrixP2 = matrixP2(:,I);
matrixDSlong2 = matrixDSlong2(:,I);
matrixDSshort2 = matrixDSshort2(:,I);
matrixDSall2 = matrixDSall2(:,I);
matrixDSsbe2 = matrixDSsbe2(:,I);

%%
fullfigure
subplot(4,1,1)
pcolor(matrixi2,matrixP2,abs(matrixDSlong2)); shading flat
axis ij
cmocean('amp')
colorbar
caxis([0 0.04])
title('Amplitude of long-term correction on RBRargo^3 CTD')
ylabel('P [dbar]')
ylim([0 100])
yyaxis right
plot(matrixi2(1,:),dpdt,'k','linewidth',2)
set(gca,'ycolor','k')
ylabel('V_p [dbar/s]')


subplot(4,1,2)
pcolor(matrixi2,matrixP2,abs(matrixDSshort2)); shading flat
axis ij
cmocean('amp')
colorbar
caxis([0 0.04])
title('Amplitude of short-term correction on RBRargo^3 CTD')
ylabel('P [dbar]')
ylim([0 100])
yyaxis right
plot(matrixi2(1,:),dpdt,'k','linewidth',2)
set(gca,'ycolor','k')
ylabel('V_p [dbar/s]')

subplot(4,1,3)
pcolor(matrixi2,matrixP2,abs(matrixDSall2)); shading flat
axis ij
cmocean('amp')
colorbar
caxis([0 0.04])
title('Amplitude of thermal inertia corrections on RBRargo^3 CTD')
ylabel('P [dbar]')
ylim([0 100])
yyaxis right
plot(matrixi2(1,:),dpdt,'k','linewidth',2)
set(gca,'ycolor','k')
ylabel('V_p [dbar/s]')

subplot(4,1,4)
pcolor(matrixi2,matrixP2,abs(matrixDSsbe2)); shading flat
axis ij
cmocean('amp')
colorbar
caxis([0 0.04])
title('Amplitude of thermal inertia corrections on SBE41CP CTD')
ylabel('P [dbar]')
xlabel('Profile number')
ylim([0 100])
yyaxis right
plot(matrixi2(1,:),dpdt,'k','linewidth',2)
set(gca,'ycolor','k')
ylabel('V_p [dbar/s]')
fsize(18)

set(gcf,'color','w')
export_fig -r200 SBE_vs_RBR_thermal_mass.png
