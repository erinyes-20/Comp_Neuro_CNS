%This module plots the frequency vs external current I and displays various thresholds of different model behavior

%Time is in msecs, voltage in mvs, conductances in m mho/mm^2, capacitance in uF/mm^2

% Accurate l1 threshold value of current is 0.0223

gkmax=.36;
vk=-77; 
gnamax=1.20;
vna=50; 
gl=0.003;
vl=-54.387; 
cm=.01; 
l1=0;           %Thresholds l1, l2, l3 are initialized
l2=0;
l3=0;

Is = 0:0.005:1;         %This is the range of current values used (resolution = dI = 0.005)

peakcount = [];         %This empty array will store frequency of spikes for every Iext in the range Is, injected into the model

for curiter=Is

    ImpCur = curiter;

    dt=0.01;
    niter=10000;
    t=(1:niter)*dt;
    iapp=ImpCur*ones(1,niter);

    v=-64.9964;
    m=0.0530;
    h=0.5960;
    n=0.3177;

    gnahist=zeros(1,niter);
    gkhist=zeros(1,niter);
    vhist=zeros(1,niter);
    mhist=zeros(1,niter);
    hhist=zeros(1,niter);
    nhist=zeros(1,niter);

    for iter=1:niter
      gna=gnamax*m^3*h; 
      gk=gkmax*n^4; 
      gtot=gna+gk+gl;
      vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
      tauv = cm/gtot;
      v=vinf+(v-vinf)*exp(-dt/tauv);
      alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
      betam = 4*exp(-0.0556*(v+65));
      alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
      betan = 0.125*exp(-(v+65)/80);
      alphah = 0.07*exp(-0.05*(v+65));
      betah = 1/(1+exp(-0.1*(v+35)));
      taum = 1/(alpham+betam);
      tauh = 1/(alphah+betah);
      taun = 1/(alphan+betan);
      minf = alpham*taum;
      hinf = alphah*tauh;
      ninf = alphan*taun;
      m=minf+(m-minf)*exp(-dt/taum);
      h=hinf+(h-hinf)*exp(-dt/tauh);
      n=ninf+(n-ninf)*exp(-dt/taun);
      vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
    end
    
    peakloc = findpeaks(vhist,'MinPeakProminence',5, 'MinPeakHeight',5);            %Finding number of spikes in 100ms
    peakcount = [peakcount, numel(peakloc)];            %Appending current peakcount to the array
    
    if numel(peakloc)==1 && l1==0           %Detecting first AP threshold l1
        l1=curiter;
    elseif numel(peakloc) > 4 && l2==0      %Detecting onset of limit cycle behavior l2
        l2=curiter;
    end

    if length(peakcount) > 4 && l3==0       %Detecting end of limit cycle and frequency decline l3
        prevpeak = peakcount(length(peakcount)-1);
        if numel(peakloc) < prevpeak
            l3=curiter;
        end
    end

end

plot(Is, peakcount)
hold on
xline(l1,'r', 'l1')
xline(l2,'r', 'l2')
xline(l3,'r', 'l3')
title('Spike frequency vs Iext curve')
xlabel('Input current')
ylabel('Spike frequency')
hold off

disp("Thus using a uniform current step of 0.005 uA, we obtain the following threshold values: ")
fprintf('First AP seen at %s uA \n',l1);
fprintf('Limit cycles initiated at %s uA \n',l2);
fprintf('Spike frequency decline seen at %s uA \n',l3);
