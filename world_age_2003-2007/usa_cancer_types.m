usa=xlsread('38408099.csv','D:E');
S=xlsread('latency&driver number.xlsx','D2:D19');
T=xlsread('latency&driver number.xlsx','C2:C19');
[~,name]=xlsread('latency&driver number.xlsx','A2:A19');
year_divide=xlsread('latency&driver number.xlsx','E2:E19');
usa_risk=usa(:,1)./usa(:,2);



for i=1:9272
    usa_Risk(i)=sum(usa_risk(floor(i/19)*19+1:i));
end
usa_Risk=usa_Risk';

t=[4 9 14 19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
t1=[14 19 24 29 34 39 44 49 54 59 64 69 74 79 84 89]';

%colon Adenocarcinoma
x=(1.2:2);
t=[19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
[p1]=polyfit(log10(t),log10(usa_Risk(4+41*38:17+41*38)),1);
plot(log10(t),log10(usa_Risk(4+41*38:17+41*38)),'o',x,polyval(p1,x),'-');
hold on
[p2]=polyfit(log10(t),log10(usa_risk(4+41*38:17+41*38)),1);
plot(log10(t),log10(usa_risk(4+41*38:17+41*38)),'o',x,polyval(p2,x),'-');

%Lung Adenocarcinoma
x=(0:2);
t=[19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
[p1]=polyfit(log10(t),log10(usa_Risk(4+78*38:17+78*38)),1);
plot(log10(t),log10(usa_Risk(4+78*38:17+78*38)),'o',x,polyval(p1,x),'-');
hold on
[p2]=polyfit(log10(t),log10(usa_risk(4+78*38:17+78*38)),1);
plot(log10(t),log10(usa_risk(4+78*38:17+78*38)),'o',x,polyval(p2,x),'-');

%Lung Adenocarcinoma
x=(9.5:10.5);
t=[19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
D=S(9)*2-2+year_divide(9).*S(9).*t;
[p]=polyfit(log10(D),log10(usa_Risk(4+78*38:17+78*38)),1);
plot(log10(D),log10(usa_Risk(4+78*38:17+78*38)),'o',x,polyval(p,x),'-');

%Acute Myeloid leukaemia
x=(10.2:11.5);
t=[19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
D=S(12)*2-2+year_divide(12).*S(12).*t;
[p]=polyfit(log10(D),log10(usa_Risk(4+232*38:17+232*38)),1);
plot(log10(D),log10(usa_Risk(4+232*38:17+232*38)),'o',x,polyval(p,x),'-');

%Acute Myeloid leukaemia
x=(1.2:2);
t=[19 24 29 34 39 44 49 54 59 64 69 74 79 84]';
[p1]=polyfit(log10(t),log10(usa_Risk(4+232*38:17+232*38)),1);
plot(log10(t),log10(usa_Risk(4+232*38:17+232*38)),'o',x,polyval(p1,x),'-');
hold on
[p2]=polyfit(log10(t),log10(usa_risk(4+232*38:17+232*38)),1);
plot(log10(t),log10(usa_risk(4+232*38:17+232*38)),'o',x,polyval(p2,x),'-');



