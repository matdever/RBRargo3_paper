clear all
close all

%% Argo CTD (32 Hz)
filename = '060668_20190815_1736.rsk';

rsk= RSKopen(filename,'readHiddenChannels',true);
rsk = RSKreaddata(rsk);

%% some post-processing
rsk = RSKcorrecthold(rsk,'channel',{'pressure','conductivity','temperature'},...
                         'action','interp');
rsk = RSKderiveseapressure(rsk);
rsk = RSKderivesalinity(rsk);
RSKprintchannels(rsk)


%% profiles
rsk = RSKtimeseries2profiles(rsk);
rsk = RSKremovecasts(rsk,'direction','up');
rsk = RSKderivedepth(rsk);
rsk = RSKderivevelocity(rsk);

% clf;RSKplotprofiles(rsk,'channel','temperature')
rsk = RSKtrim(rsk,'reference','sea pressure','range',[-inf 1],'action','remove');
rsk = RSKtrim(rsk,'reference','sea pressure','range',[9.5 inf],'action','remove');

vcol = getchannelindex(rsk,'velocity');
pcol = getchannelindex(rsk,'sea pressure');
for k=1:length(rsk.data),
    
    jj = rsk.data(k).values(:,pcol) > 2 & ...
         rsk.data(k).values(:,pcol) < 4;
    vave = 100*nanmean(rsk.data(k).values(jj,vcol));
    
    %rsk.data(k).descentrate = sprintf('%2.1f',vave);    
    rsk.data(k).descentrate = vave;
    
    rsk.data(k).mtime = rsk.data(k).tstamp(1);
    rsk.data(k).time = datestr(rsk.data(k).tstamp(1),'HH:MM:SS');

end

[~,ind] = sort([rsk.data.descentrate]);
rsk.data = rsk.data(ind);



%% apply an advance to temperature
deltat = -0.35;
nsamp = round(deltat/readsamplingperiod(rsk));
crsk = RSKalignchannel(rsk,'channel','Temperature','lag',nsamp);
crsk = RSKderivesalinity(crsk);

%% sharpen thermistor with Fozdar
crsk2 = crsk;
crsk2 = RSKcorrecttau(crsk,'channel','temperature','tauresponse',0.35);
crsk2 = RSKderivesalinity(crsk2);

%% plot profiles

fullfigure
chans = {'temperature','conductivity','salinity'};
[h1,ax] = RSKplotprofiles(crsk,'channel',chans,'profile',6);
arrayfun(@(x) hold(x,'on'),ax)
[h2,~] = RSKplotprofiles(rsk,'channel',chans,'profile',6);
%[h3,~] = RSKplotprofiles(crsk2,'channel',chans,'profile',6);

set(h1,'linestyle','-','color','r','linewidth',1.5)
set(h2,'linestyle','-','color','k','linewidth',1.5)
ylim([1 6])

legend([h2(1) h1(1)],'raw','corrected for C-T lag')
fsize(30)

set(ax(1),'xlim',[20 35])
set(ax(3),'xlim',[0 17])

axes(ax(1)); xlabel('Temperature [°C]'); title(''); ylabel('Pressure [dbar]')
axes(ax(2)); xlabel('Conductivity [mS/cm]'); title(''); ylabel('')
axes(ax(3)); xlabel('Salinity [PSU]'); title(''); ylabel('')

set(ax,'xaxislocation','top')

set(gcf,'color','w')
export_fig('Fig8_Profile_example_CTlag.png','-r200')