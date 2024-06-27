
pC=0.3     //Assortment probability for cooperators
pN=0.33    //Assortment probability for noncooperators
LR=0.2     //Intensity of leaving the receiver state
LD=0.4     //Intensity of leaving the donor state
B=0.3;      //Benefit
C=0.1;      //Cost
D=1-B-C;
S=B+C+D;

q=[0:0.02:1];
p=[0:0.02:1];


  qC=(LR+LD-(1-p).*B-D+sqrt((LR+LD-(1-p).*B-D)^2+4.*LR.*((1-p).*B+D)))./(2.*((1-p).*B+D))

qCpC=(LR+LD-(1-pC)*B-D+sqrt((LR+LD-(1-pC)*B-D)^2+4*LR.*((1-pC).*B+D)))./(2*((1-pC)*B+D))

  qN=(LR+LD-(1-p).*B-C-D+sqrt((LR+LD-(1-p)*B-C-D)^2+4*LR*(B.*(1-p)+C+D)))./(2*((1-p)*B+C+D))

qNpN=(LR+LD-(1-pN)*B-C-D+sqrt((LR+LD-(1-pN)*B-C-D)^2+4*LR*(B.*(1-pN)+C+D)))./(2*((1-pN)*B+C+D))


  //     (1-(LR+LD-(1-p).*B-D+sqrt((LR+LD-(1-p).*B-D)^2+4.*LR.*((1-p).*B+D)))./(2.*((1-p).*B+D))).*((1-p).*B+D)+C
mortqC=(1-((LR+LD-(1-p).*B-D+sqrt((LR+LD-(1-p).*B-D)^2+4.*LR.*[(1-p).*B+D]))./(2.*[(1-p).*B+D]))).*((1-p).*B+D)+C
mortqD=(1-((LR+LD-(1-p).*B-C-D+sqrt((LR+LD-(1-p).*B-C-D)^2+4.*LR.*(B.*(1-p)+C+D)))./(2.*((1-p).*B+C+D)))).*((1-p).*B+D+C)

if qCpC>1 then qCpC=1
end
if qCpC<0 then qCpC=0
end
if qNpN>1 then qNpN=1
end
if qNpN<0 then qNpN=0
end



deff('f=fC(q,p)','f=(1-q)*((1-p)*B+D)+C')
deff('f=fN(q,p)','f=(1-q)*((1-p)*B+D+C)')
//deff('r=f(qn,qc)','r=(1-qn-qc)./(2-3.*qn-qc+0.0000000000000000000001)');

fitC=fC(qCpC,pC)
fitN=fN(qNpN,pN)

fitdif=fitC-fitN


f0=scf(1)
clf
f0.color_map=graycolormap(8) ;

fplot3d1(q,p,fC,flag=[1,3,4],leg="q@p@mortality",ebox=[0,1,0,1,0,S])
fplot3d1(q,p,fN,flag=[1,3,4],leg="q@p@mortality",ebox=[0,1,0,1,0,S])

f1=scf(2)
clf
f1.color_map=graycolormap(8) ;

fplot3d1(q,p,fC,flag=[1,3,4],leg="q@p@mortality",ebox=[0,1,0,1,0,S])

f2=scf(3)
clf
f2.color_map=graycolormap(8) ;

fplot3d1(q,p,fN,flag=[1,3,4],leg="q@p@mortality",ebox=[0,1,0,1,0,S])

f3=scf(4)
clf
f3.color_map=graycolormap(8) ;

plot2d(q,fC(q,pC), style=1)
//plot2d(q,fC(q,pN),style=1)
plot2d(q,fN(q,pN), style=7)
//plot2d(q,fN(q,pC),style=4)
plot2d(p*0+qCpC , q , style=1, rect=[0,0,1,1-B])
plot2d(p*0+qNpN , q , style=7, rect=[0,0,1,1-B])
xlabel("q Donor", "fontsize", 2);
ylabel("mortality", "fontsize", 2);
legends(['Cooperators';'Noncooperators'],[1,7],opt="?")
xtitle('relative mortality');

f4=scf(5)
clf
f4.color_map=graycolormap(8) ;

plot2d(p,qC, style=1 ) //, rect=[0,0,1,1])
//plot2d(q,fC(q,pN),style=1)
plot2d(p,qN, style=7 ) //, rect=[0,0,1,1])
xlabel("assortment probability", "fontsize", 2);
ylabel("stable fraction of Donors", "fontsize", 2);
legends(['Cooperators';'Noncooperators'],[1,7],opt="?")
xtitle('plot of stable Donor/receiver role distributions described as the fraction of Donors');

f5=scf(6)
clf
f5.color_map=graycolormap(8) ;

plot2d(p,mortqC, style=1)
//plot2d(q,fC(q,pN),style=1)
plot2d(p,mortqD, style=7)
//plot2d(q,fN(q,pC),style=4)
plot2d(p*0+pC , q , style=1, rect=[0,0,1,1-B])
plot2d(p*0+pN , q , style=7, rect=[0,0,1,1-B])
xlabel("assortment probability", "fontsize", 2);
ylabel("mortality", "fontsize", 2);
legends(['Cooperators';'Noncooperators'],[1,7],opt="?")
xtitle('relative mortality with substituted q');





