clear all; %close all
% dependencies:
% 1 - GSW
% 2 - SBE-toolbox
% 3 - RSK-tools

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates colors for RBR and SBE
rbrred = [186 12 47]/256;
sbeblue = [17 25 92]/256;

% plot both results for default and custom pressure correction ceofs.
for custom = [0,1]
    
    
    %% SBE data
    % Import data
    filename = 'datasets/in2019_v06_004.cnv';
    sbe = readSBScnv(filename);
    
    % compute Salinity and Conservative temperature
    SBESP = gsw_SP_from_C(sbe.c0Sm*10,sbe.t090C,sbe.pm);
    [SBESA, ~] = gsw_SA_from_SP(SBESP,sbe.pm,120,-14);
    SBECT = gsw_CT_from_t(SBESA,sbe.t090C,sbe.pm);
    SBEp = sbe.pm-0.44;
    SBE.time = sbe.timeS/86400+datenum(sbe.instrumentheaders.NMEAUTCTime);
    
    SBESA = SBESA(SBE.time<datenum('22-Oct-2019 12:59:35'));
    SBECT = SBECT(SBE.time<datenum('22-Oct-2019 12:59:35'));
    SBEp = SBEp(SBE.time<datenum('22-Oct-2019 12:59:35'));
    SBE.time = SBE.time(SBE.time<datenum('22-Oct-2019 12:59:35'));
    
    % Bin the data in pressure (bin size = Prez)
    Prez = 10;
    Pgrid = 0 : Prez : 2000;
    for pp = 1:length(Pgrid)
        try
            SBE.CTgrid(pp) = nanmean(SBECT(SBEp>Pgrid(pp)-Prez/2 & SBEp<Pgrid(pp)+Prez/2));
            SBE.SAgrid(pp) = nanmean(SBESA(SBEp>Pgrid(pp)-Prez/2 & SBEp<Pgrid(pp)+Prez/2));
        catch
            SBE.CTgrid(pp) = NaN;
            SBE.SAgrid(pp) = NaN;
        end
    end
    
    % Bin the data in temperature (bin size = Trez)
    Trez = 0.05;
    Tgrid = min(SBE.CTgrid)+Trez/2:Trez:max(SBE.CTgrid)-Trez/2;
    for pp = 1:length(Tgrid)
        SBE.SAgrid2(pp,1) = nanmean(SBE.SAgrid(SBE.CTgrid>Tgrid(pp)-Trez/2 & SBE.CTgrid<Tgrid(pp)+Trez/2));
    end; clear pp
    
    
    %% PLot SBE data 
    % Create the T-S plot with SBE data
    if custom == 0
        f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if custom
        subplot(5,3,[2 5 8])
    else
        subplot(5,3,[1 4 7])
    end
    h1 = plot(SBE.SAgrid,SBE.CTgrid,'b','linewidth',2);
    hold on
    % Adds point to place depth levels on T-S plot for reference
    for zz = 800:200:2000
        [~,I] = min(abs(Pgrid-zz));
        scatter(SBE.SAgrid(I),SBE.CTgrid(I),100,'ok','filled')
        text(SBE.SAgrid(I)-0.01,SBE.CTgrid(I),[num2str(zz),' m'],'fontsize',16,'horizontalalignment','right')
    end; clear zz
    clearvars -except rbrred sbeblue Prez Pgrid h1 SBE f1 custom
    
    
    %% RBR
    % Load RBR file
    for file = 1:2
        if file == 1
            filename = 'datasets/060669_20191026_0201.rsk';
        elseif file == 2
            filename = 'datasets/060671_20191026_0223.rsk';
        end
        
        % Open RSK
        rsk = RSKopen(filename,'readHiddenChannels',true);
        rsk = RSKreaddata(rsk);
        Fs = round(1./rsk.continuous.samplingPeriod*1000);
        rsk = RSKreadcalibrations(rsk);
        RSKprintchannels(rsk)
        
        % Update pressure correction from linear (depreciated) to cubic if
        % necessary
        if rsk.calibrations(1).x2 == 4e-7
            warning('Updating coeficients from linear')
            rsk = update_pressure_correction_lin2cubic(rsk);
        elseif rsk.calibrations(1).x2 == 1.8732e-06
            warning('Updating coeficients')
            if custom
                rsk = update_pressure_correction_default2custom(rsk);
            else
                myp = rsk.data.values(:,3) -  rsk.calibrations(1).x6;
                rsk.data.values(:,1) = rsk.data.values(:,1).*(1 + rsk.calibrations(1).x2.*myp + rsk.calibrations(1).x3.*myp.^2 + rsk.calibrations(1).x4.*myp.^3);
            end
            
        end
        
        % Select a profile for comparison
        for ii = 4
            
            % Crop data to that downcast
            ind = find(rsk.data.tstamp>=rsk.profiles.downcast.tstart(ii) & rsk.data.tstamp<=rsk.profiles.downcast.tend(ii));
            time = rsk.data.tstamp(ind);
            if file == 2
                % Implement the K-factor derived by comparing with bottle
                % samples.
                C = rsk.data.values(ind,1)*0.9997;
            elseif file == 1
                % Implement the K-factor derived by comparing with bottle
                % samples.
                C = rsk.data.values(ind,1)*1.0001;
            else
                C = rsk.data.values(ind,1);
            end
            t = rsk.data.values(ind,2);
            
            idx = getchannelindex(rsk,'CT Cell Temperature');
            Tcond = rsk.data.values(ind,idx);
            clear idx
            
            % compute sea pressure
            if file == 1
                p = rsk.data.values(ind,3)-10.11;
            elseif file == 2
                p = rsk.data.values(ind,3)-10.16;
            end
            
            % compute Salinity and Conservative temperature
            SP = gsw_SP_from_C(C,t,p);
            [SA, ~] = gsw_SA_from_SP(SP,p,120,-14);
            CT = gsw_CT_from_t(SA,t,p);
            
            % Bin the data on pressure bins (bin size = Prez)
            for pp = 1:length(Pgrid)
                try
                    CTgrid(pp) = nanmean(CT(p>Pgrid(pp)-Prez/2 & p<Pgrid(pp)+Prez/2));
                    SAgrid(pp) = nanmean(SA(p>Pgrid(pp)-Prez/2 & p<Pgrid(pp)+Prez/2));
                catch
                    CTgrid(pp) = NaN;
                    SAgrid(pp) = NaN;
                end
            end
            
            % Bin the data on temperature bins (bin size = Trez)
            Trez = 0.05;
            Tgrid = min(SBE.CTgrid)+Trez/2:Trez:max(SBE.CTgrid)-Trez/2;
            for pp = 1:length(Tgrid)
                SAgrid2(pp,1) = nanmean(SAgrid(CTgrid>Tgrid(pp)-Trez/2 & CTgrid<Tgrid(pp)+Trez/2));
            end; clear pp
            %%
            
            figure(f1)
            if custom
                subplot(5,3,[2 5 8])
            else
                subplot(5,3,[1 4 7])
            end
            if file == 1
                h2(file) = plot(SAgrid,CTgrid,'color',rbrred,'linestyle','-','linewidth',2);
            elseif file == 2
                h2(file) = plot(SAgrid,CTgrid,'color',rbrred,'linestyle','--','linewidth',2);
            end
            
            if file == 1
                %f2 = figure('Position',[869 398 877 360]);
                if custom
                    subplot(5,3,[11 14])                
                else
                    subplot(5,3,[10 13])
                end
                if custom == 0
                    h3(file) = histogram(SAgrid(Pgrid>800)-SBE.SAgrid(Pgrid>800),-0.06:0.0005:0.06,'facecolor',rbrred,'facealpha',.8);
                else
                    h3(file) = histogram(SAgrid(Pgrid>800)-SBE.SAgrid(Pgrid>800),-0.06:0.0005:0.06,'facecolor',rbrred,'facealpha',.8);
                end
                %h3(file) = histogram(SBE.SAgrid2(Tgrid<SBE.CTgrid(Pgrid==800)) - SAgrid2(Tgrid<CTgrid(Pgrid==800)),-0.036:0.0005:0.036,'facecolor',rbrred,'facealpha',.8);
                hold on
            elseif file == 2
                if custom
                    subplot(5,3,[11 14])
                else
                    subplot(5,3,[10 13])
                end
                if custom == 0
                    h3(file) = histogram(SAgrid(Pgrid>800)-SBE.SAgrid(Pgrid>800),-0.06:0.0005:0.06,'facecolor',rbrred,'facealpha',.5);
                else
                    h3(file) = histogram(SAgrid(Pgrid>800)-SBE.SAgrid(Pgrid>800),-0.06:0.0005:0.06,'facecolor',rbrred,'facealpha',.5);
                end
                %h3(file) = histogram(SBE.SAgrid2(Tgrid<SBE.CTgrid(Pgrid==800)) - SAgrid2(Tgrid<CTgrid(Pgrid==800)),-0.036:0.0005:0.036,'facecolor',rbrred,'facealpha',.5);
            end
            
        end
    end
    
    figure(f1)
    if custom
        subplot(5,3,[2 5 8])
        text(34.605,8.36,'C.','fontsize',18,'fontweight','b')
    else
        subplot(5,3,[1 4 7])
        text(34.605,8.36,'A.','fontsize',18,'fontweight','b')
    end
    % Make plot pretty
    xlim([34.6 35.1])
    ylim([2 8.2])
    [X,Y] = meshgrid(32:0.1:37,1:0.1:20);
    rho = gsw_rho(X,Y,0);
    [C,h] = contour(X,Y,rho,1020:.2:1040,'k');
    clabel(C,h)
    xlabel('Absolute salinity [g/kg]')
    ylabel('Conservative temperature [°C]')
    legend([h1 h2],'SBE9','SN060669','SN060671','SN060672','location','northwest')
    set(gca,'fontsize',16)
    axis square
    
    %title('YMC (2019)')
    set(gcf,'color','w')
    
    if custom
        subplot(5,3,[11 14])
        text(-0.064,36.5,'D.','fontsize',18,'fontweight','b')        
    else
        subplot(5,3,[10 13])
        text(-0.064,36.5,'B.','fontsize',18,'fontweight','b')                
    end
    % Make plot pretty
    xlabel('Salinity difference (RBR - SBE) [g/kg]')
    ylabel('PDF')
    set(gca,'fontsize',16)
    %title('YMC (2019)')
    set(gcf,'color','w')
    grid on; grid minor
    ylim([0 35])
    YLIM = get(gca,'ylim');
    patch([-0.006 -0.006 0.006 0.006 -0.006],[YLIM(1) YLIM(2) YLIM(2) YLIM(1) YLIM(1)],[.7 .7 .7],'facealpha',.3)
    legend(h3,'SN060669','SN060671','SN060672','location','northwest')
    %mine = get(gca,'position');
    %mine(4) = 0.22;
    %set(gca,'position',mine)
    set(gcf,'color','w')

    
end

export_fig -r200 compression_errors_YMC.png

%% FUNCTION update_pressure_correction
function rsk_corrected = update_pressure_correction_lin2cubic(rsk)


rsk_corrected = rsk;


% linear coefficient for pressure correction
x2linear = rsk.calibrations(1).x2;

% calibration absolute pressure
x4linear = rsk.calibrations(1).x4;

% conductivity correction coefficients
if rsk.instruments.serialID == 60669
    disp('Coefs for 060669')
    x2 =  1.5e-06;
    x3 = -5.28e-10;
    x4 =  7.76e-14;
elseif rsk.instruments.serialID == 60671
    disp('Coefs for 060671')
    x2 =  1.57e-06;
    x3 = -4.42e-10;
    x4 =  5.37e-14;
else
    x2 =  1.8732e-06;
    x3 = -7.7689e-10;
    x4 =  1.4890e-13;
end

fp = @(P) x2.*P + x3.*P.^2 + x4.*P.^3;


% find indices to extract pressure and conductivity from data table
pcol = getchannelindex(rsk,'pressure');
ccol = getchannelindex(rsk,'conductivity');


% loop over profiles and compute properly corrected conductivity
for k=1:length(rsk.data)
    
    c_orig = rsk.data(k).values(:,ccol);
    p = rsk.data(k).values(:,pcol) - x4linear;
    
    c_uncorrected = c_orig.*(1 + x2linear.*p);
    
    c_corrected = c_uncorrected./(1 + fp(p));
    
    rsk_corrected.data(k).values(:,ccol) = c_corrected;
    
end


% udpate calibration table with new coefficients
rsk_corrected.calibrations(1).x2 = 'cubic';
end


function rsk_corrected = update_pressure_correction_default2custom(rsk)


rsk_corrected = rsk;


% linear coefficient for pressure correction
x2default = rsk.calibrations(1).x2;
x3default = rsk.calibrations(1).x3;
x4default = rsk.calibrations(1).x4;

% calibration absolute pressure
x6default = rsk.calibrations(1).x6;

% conductivity correction coefficients
if rsk.instruments.serialID == 60669
    disp('Coefs for 060669')
    x2 =  1.5e-06;
    x3 = -5.28e-10;
    x4 =  7.76e-14;
elseif rsk.instruments.serialID == 60671
    disp('Coefs for 060671')
    x2 =  1.41e-06;
    x3 = -5.25e-10;
    x4 =  7.90e-14;
else
    x2 =  1.8732e-06;
    x3 = -7.7689e-10;
    x4 =  1.4890e-13;
end

fp = @(P) x2.*P + x3.*P.^2 + x4.*P.^3;


% find indices to extract pressure and conductivity from data table
pcol = getchannelindex(rsk,'pressure');
ccol = getchannelindex(rsk,'conductivity');


% loop over profiles and compute properly corrected conductivity
for k=1:length(rsk.data)
    
    c_orig = rsk.data(k).values(:,ccol);
    p = rsk.data(k).values(:,pcol) - x6default;
    
    c_uncorrected = c_orig.*(1 + x2default.*p + x3default.*p.^2 + x4default.*p.^3);
    
    c_corrected = c_uncorrected./(1 + fp(p));
    
    rsk_corrected.data(k).values(:,ccol) = c_corrected;
    
end


% udpate calibration table with new coefficients
%rsk_corrected.calibrations(1).x2 = 'cubic';
end
