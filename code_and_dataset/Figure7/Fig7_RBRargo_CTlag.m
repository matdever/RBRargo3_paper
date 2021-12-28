clear
% This code computed the optimal C-T lag using all profiles collected by 6
% RBRConcerto3 mounted on rosettes (downcasts only), during 5 different
% cruises

profnum = 0;
for file = 1:20
    % Loop through all the files names to be analyzed
    if file == 1
        filename = 'datasets/060667_20191122_0900.rsk';
    elseif file == 2
        filename = 'datasets/060667_20181126_1104.rsk';
        
    elseif file == 3
        filename = 'datasets/060668_20191122_0852.rsk';
    elseif file == 4
        filename = 'datasets/060668_20181126_1107.rsk';
        
    elseif file == 5
        filename = 'datasets/060669_20191026_0201.rsk';
    elseif file == 6
        filename = 'datasets/060669_20191028_1241.rsk';
    elseif file == 7
        filename = 'datasets/060669_20191030_0005.rsk';
    elseif file == 8
        filename = 'datasets/060669_20191103_1237.rsk';
    elseif file == 9
        filename = 'datasets/060669_20191104_1718.rsk';
    elseif file == 10
        filename = 'datasets/060670_20191122_0907.rsk';
        
    elseif file == 11
        filename = 'datasets/060671_20181126_1102.rsk';
    elseif file == 12
        filename = 'datasets/060671_20191026_0223.rsk';
    elseif file == 13
        filename = 'datasets/060671_20191028_1237.rsk';
    elseif file == 14
        filename = 'datasets/060671_20191030_0015.rsk';
    elseif file == 15
        filename = 'datasets/060671_20191103_1229.rsk';
    elseif file == 16
        filename = 'datasets/060671_20191104_1720.rsk';
    elseif file == 17
        filename = 'datasets/060671_20191105_1959.rsk';
    elseif file == 18
        filename = 'datasets/060669_20180808_0329.rsk';
    elseif file == 19
        filename = 'datasets/060671_20180808_0327.rsk';
    elseif file == 20
        filename = 'datasets/060672_20180614_2233.rsk';
    end
    file
    
    % Open RSK-file and load profiles, just in case they are not already
    % identified in the file
    rsk = RSKopen(filename);
    rsk = RSKreaddata(rsk);
    rsk = RSKfindprofiles(rsk);
    
    % Extract sampling frequency [in Hz]
    Fs = round(1./rsk.continuous.samplingPeriod*1000);
    
    % start counters to keep track of the number of segments and how many
    % are ignored
    if exist('segnum','var')==0
        segnum = 1;
        skipped = 0;
    end
    
    % for each downcast
    for ii = 1:length(rsk.profiles.downcast.tstart)
        
        % Extract the data of the profile
        ind = find(rsk.data.tstamp>=rsk.profiles.downcast.tstart(ii) & rsk.data.tstamp<=rsk.profiles.downcast.tend(ii));
        time = rsk.data.tstamp(ind);
        C = rsk.data.values(ind,1);
        T = rsk.data.values(ind,2);
        P = rsk.data.values(ind,3);
        
        % compute difference
        Cdiff = diff(C);
        Tdiff = diff(T);
        
        % Exclude profiles that are too shallow or strange (e.g., don't start
        % at the surface)
        if min(P)>50
            disp(['Skip cast ',num2str(ii),' - weird profile'])
            continue
        end
        
        profnum = profnum + 1;
        % 1st segment is at 50 m deep to avoid ML where lag will be skewed
        % towards 0
        Pstart = find(P>50,1,'first');
        Pend = find(smooth(diff(P(Pstart:end))./diff(time(Pstart:end)*86400),30*Fs)>0.25,1,'last');
        Pend = Pend + Pstart - 1 ;
        segment = 5*Fs; % 5s segments
        
        % For each 5s segment with 50% overlap
        for tt = Pstart:segment/2:Pend-segment
            [thecov,lags] = xcov(Tdiff(tt:tt+segment),Cdiff(tt:tt+segment),50,'coef');
            
            % If the covariance is too weak, segment is rejected, Nans are
            % recorded
            if max(thecov)<0.5
                warning(['Cov is weak (',num2str(max(thecov)),')'])
                lag(segnum) = NaN;
                lag_s(segnum) = NaN;
                dpdt(segnum) = NaN;
                segnum = segnum+1;
                skipped= skipped + 1;
                continue
            else
                % take maximum covariance
                [C,I] = max(thecov);
                if I == 1 || I == length(thecov)
                    warning(['Cov is weak (',num2str(max(thecov)),')'])
                    lag(segnum) = NaN;
                    dpdt(segnum) = NaN;
                    segnum = segnum+1;
                    skipped = skipped + 1;
                    continue
                end
                % fit 2nd order polynomial
                f = fit(lags(I-1:I+1)',thecov(I-1:I+1),'poly2');
                % Get maximum of the fit
                lag(segnum) = -f.p2./(2*f.p1);
                lag_s(segnum) = lag(segnum)./Fs;
                % record corresponding fall-rate
                dpdt(segnum) = nanmean(diff(P(tt:tt+segment))./diff(time(tt:tt+segment)*86400));
                
            end
            segnum = segnum+1;
        end
    end
end

%%
% display the total number of segments used.
disp(['N = ',num2str(segnum - skipped)])
N1 = segnum - skipped;
figure
h5 = histogram(lag_s(lag_s>0.15),0.16:0.01:1,'facecolor','k','Normalization','probability','facealpha',.3);
hold on
% Add Gaussian fit
[myN,myE] = histcounts(lag_s(lag_s>0.15 & lag_s<.6),0.16:0.01:1,'Normalization','probability');
[mdl,gof] = fit(((myE(1:end-1)+myE(2:end))/2)',myN','gauss1');

%mdl = createFit(lag_s(lag_s>0.15 & lag_s<.55));
h6 = plot(0.16:0.01:1,mdl.a1*exp(-(([0.16:0.01:1] - mdl.b1)./mdl.c1).^2),'-r','linewidth',2)
%plot(0.16:0.01:1,normpdf(0.16:0.01:1,mdl.mu,mdl.sigma),'-r')
h7 = line([mdl.b1 mdl.b1], get(gca,'ylim'),'color','k','linewidth',2);
%set(h6,'color','k','linewidth',2)
xlabel('Lag [s]')
ylabel('PDF')
%h7 = line([0.7 0.7],get(gca,'ylim'),'color','r','linewidth',2);
%h8 = line([1 1],get(gca,'ylim'),'color','r','linewidth',2,'linestyle','--');
legend([h6 h7],...
    ['Gaussian fit (R^2 = ',num2str(gof.rsquare,'%0.2f'),'; N = ',num2str(N1),')'],...
    ['\mu = ',num2str(mdl.b1,'%0.2f'),' s; \sigma =  ',num2str(mdl.c1*sqrt(2),'%0.2f'),' s'],...
    'location','northoutside')
set(gca,'fontsize',16)
grid on; grid minor
xlim([0.15 .6])

set(gcf,'color','w')
export_fig -r200 RBRargo_CTlag_PDF.png
%%
figure
histogram2(dpdt(lag_s>0.15 & lag_s<0.6),lag_s(lag_s>0.15 & lag_s<0.6),0:.05:2,0.1:0.01:0.7,'FaceColor','flat')
ylim([0.15 .6])
view(0,90)
xlabel('Fall rate [dbar/s]')
ylabel('Lag [s]')
set(gca,'fontsize',16)
%title({'PDF of lags as a function of fall rate for slow thermistors','(SN060667, SN060668, SN060669, SN060670)'})
hold on

% Add average of guassian fit
%mdl = createFit(lag_s(lag_s>0.15 & lag_s<0.6));
h1 = line(get(gca,'xlim'),[mdl.b1 mdl.b1],[max(get(gca,'zlim')) max(get(gca,'zlim'))],'linewidth',3,'color','r');

% Add 25, 50, and 75th percentiles
[N,X,Y] = histcounts2(dpdt(lag_s>0.15 & lag_s<0.6),lag_s(lag_s>0.15 & lag_s<0.6),0:.05:2,0.1:0.01:0.7);
for nn = 1:size(N,1)
    ind = find(lag_s>0.15 & lag_s<0.6 & dpdt>=X(nn) & dpdt<X(nn+1));
    distri(nn,1) = prctile(lag_s(ind),25);
    distri(nn,2) = prctile(lag_s(ind),50);
    distri(nn,3) = prctile(lag_s(ind),75);
end
x = (X(1:end-1)+X(2:end))/2;
h2 = plot3(x,distri(:,1),150*ones(size(x)),'-k','linewidth',2);
h3 = plot3(x,distri(:,2),150*ones(size(x)),'--k','linewidth',3);
h4 = plot3(x,distri(:,3),150*ones(size(x)),'-k','linewidth',2);

legend([h1 h2 h3 h4],['Constant lag = ',num2str(mdl.b1,'%0.2f'),'s'],...
    '25th percentile',...
    'median',...
    '75th percentile')
set(gca,'color','w','zscale','log','ColorScale','log')
colorbar

text(1.5,0.17,150,['N = ',num2str(N1)],'fontsize',18,'backgroundcolor','w')

set(gcf,'color','w')
export_fig -r200 RBRargo_CTlag_2DPDF.png