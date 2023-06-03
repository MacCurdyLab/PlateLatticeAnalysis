function rawData = getRawData(fName,sName,dropHeight)

[AllData] = Extractor(fName,'matchname',sName,'dropHeight',dropHeight);

Data = {};
Data2 = {};
for i = 1:length(AllData)
    Data = [Data [AllData{i}.timeF AllData{i}.Force(:,2)]];
    Data2= [Data2 [AllData{i}.timeX AllData{i}.X]];
end

%now average these responses
[t,f] = shadedErrorPrepCell(Data);
[t2,x2] = shadedErrorPrepCell(Data2);
tRaw = linspace(0,0.05,1000);
fRaw = interp1(t,nanmean(f),tRaw);
dRaw = [tRaw' fRaw'];
m_Impact = AllData{1}.mass;
v_Impact = abs(AllData{1}.vin)*1000;
E_Impact = 0.5*m_Impact*v_Impact^2;

rawData.m = m_Impact;
rawData.v = v_Impact;

DataFX = {};
for k=1:length(AllData)
    x = AllData{k}.XSync';
    F2 = AllData{k}.FSync2';
    F1 = AllData{k}.FSync1';

    %retrigger Force / Disp curve
    Ftrigger = 5; %N
    Si = find(F1>Ftrigger,1);
    x = x(Si:end) - x(Si);
    F2 = F2(Si:end) - F2(Si);
    F1 = F1(Si:end) - F1(Si);

    DataFX = [DataFX [-x F1]];
end

[x,fx] = shadedErrorPrepCell(DataFX);

xRaw = linspace(0,0.016,1000);
fRaw = interp1(x,nanmean(fx),xRaw);
fxRaw = [xRaw' fRaw'];

rawData.fx = fxRaw;
rawData.xt = dRaw;

rawData.shadedT = t;
rawData.shadedFt = f;
rawData.shadedTX = t2;
rawData.shadedXX = x2;
rawData.shadedX = x;
rawData.shadedFx = fx;

end


