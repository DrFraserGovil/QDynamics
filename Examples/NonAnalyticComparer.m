set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',23)

%% options
global errFloor legendtile folder width axs mode minT maxT;
errFloor = -2;
axs = 22;
legendtile = 3;
width = 4;
minT = 1e-1 ;
maxT = 30;
folder = "Output/Harmonic";
baseLineFile = "Sym_Euler_2_N52.dat";
triggerList = ["Mag_Euler_0","Sym_Euler_1","Sym_Euler_2","Brute"];
mode = 'log';
%% some persistent global quantities
global labelList triggerId cols positionInterp energyInterp momentumInterp saveLines
labelList = ["Analytical Solution"];
triggerId = zeros(size(triggerList))-1;
cols = [colororder; rand(14,3)];
saveLines = {};
%% file checks
fileList = organiseFiles(folder,baseLineFile);

%% main plotting loop
preparePlot();
plotBase(baseLineFile)

for file = fileList
    disp(file);
    plotFile(file,triggerList)
end

plotSavedLines;
finalStyling();

function y =clarifyer(x)
    global errFloor;
    errorMag = 20000;
    
    err = 10^(errFloor-1);

    y = (x + err)*(1.0+rand/errorMag);
end

function envelopePlot(x,y,color,width,style,active,tile)
    global errFloor minT saveLines
    y(isnan(y)|isinf(y)) = errFloor -1;
 
    if active
%         [~,wp] = findpeaks(y);
%         [~,wt] = findpeaks(-y);
%         xp = x(wp);
%         yp = y(wp);
%         xt = x(wt);
%         yt = y(wt);
%         plot(xp,yp,'Color',color,"LineWidth",width,'LineStyle',style);
%         plot(xt,yt,'Color',color,"LineWidth",width,"HandleVisibility","off",'LineStyle',style);
          subselect = (x >= minT);
          x = x(subselect);
          y = y(subselect);
%            s = smooth(y,length(y)/100);
          s = y;
           plot(x,s,'Color',[color,0.15],"LineWidth",width/2,"HandleVisibility","off",'LineStyle',style);
%           plot(x,cummax(s),'Color',color,"LineWidth",width,'LineStyle',style);
           N = 10;

          capsule = {tile,smooth(x,N),smooth(cummax(s),N),color,width,style};
          saveLines{end+1} = capsule;
          
    else
        plot(x,y,'Color',color,"LineWidth",width,'LineStyle',style);
    end
end
function smoothPlot(x,diff,color,width,style)
    global minT
    y = log10(clarifyer(diff));
     if (x(1) < minT*0.9)
        x(1) = minT*0.9;
    end
    smoothing = 10;
    goodPoints = sum(~isnan(y));
    
    while goodPoints < smoothing
        smoothing = round(smoothing/5);
    end
    smoothing = max(1,smoothing);
    smoothed = smooth(y,smoothing);
    
    plot(x(1:end-smoothing),smoothed(1:end-smoothing),'Color',color,"LineWidth",width,'LineStyle',style);
   
end
function [r] = organiseFiles(folder,target)
    q = dir(folder);
    r = convertCharsToStrings({q.name});
    r(r == "." | r == "..") = [];
    r(r == target) = [];
end
function preparePlot()  
    f=figure(1);
    f.Units = 'normalized';
    size = [0,0,0.23,1.5];
%     size = [0,0,0.3,1];
    f.OuterPosition = size;
    f.Position = size;
%     f.Position =
    clf;
    tiledlayout(3,1,"TileSpacing","Compact","Padding","Compact");
end
function plotBase(file)
     global legendtile folder positionInterp energyInterp momentumInterp errFloor;
     nexttile(1);
    baseW = 4;
    f = readtable(folder + "/" +file);	
    hold on;
    envelopePlot(f.t,f.q0,'k',baseW,'-',false);
    hold off;
    
%     
%     nexttile(4);
%     hold on;
%     yy = log10(abs(f.H - f.H(1)) + 10^(errFloor - 1));
%     envelopePlot(f.t,yy,[0.1,0.1,0.1],baseW,'-',true);
%     hold off;
    if legendtile ~= 1
        nexttile(legendtile);
        hold on;
        plot(NaN,NaN,'Color','k','LineWidth',baseW);
        hold off;
    end
    %% prepare interpolants
    positionInterp = griddedInterpolant(f.t,f.q0);
    energyInterp = griddedInterpolant(f.t,f.H);
%     momentumInterp = griddedInterpolant(f.t,f.L);
end

function plotFile(fileName,fileList)
    global width triggerId folder cols labelList width positionInterp energyInterp momentumInterp
    typeList = ["-","--",":"];
    trig = -1;
    for i = 1:length(fileList)
        if contains(fileName,fileList(i))
            trig = i;
            triggerId(i) = triggerId(i) + 1;
        end
    end
    
	if trig > -1
        if trig > 0
            color = cols(trig,:);
            sID = mod(triggerId(trig),length(typeList))+1;
            loopies = floor(triggerId(trig)/length(typeList));
            color = color * 0.7^loopies;
            style = typeList(sID);
        else
            color = 'k';
            style = '-';
        end

        
		f = readtable(folder + "/" +fileName);	
        fprintf("\tFile Loaded....");
        if height(f) > 0
     
            labelList(end+1) = extractName(fileName);
            theta = f.q0;
            theta(isnan(theta)) = 1e-15;

            nexttile(1);
            hold on;
            width = width/2;
            envelopePlot(f.t,theta,color,width,style,false,1);
            hold off;
            width = width * 2;
            fprintf("Plot 1 Complete...");

            nexttile(2);
            hold on;
            orig = positionInterp(f.t);
            diff = abs( (theta - orig));
            envelopePlot(f.t,diff,color,width,style,true,2);
            hold off;
            fprintf("Plot 2 Complete...");

            nexttile(3);
            hold on;
            orig = f.H(1);%energyInterp(f.t);
            diff = abs((f.H - orig)./(orig + 1e-10));
            envelopePlot(f.t,diff,color,width,style,true,3);
            hold off;
            fprintf("Plot 3 Complete\n\n");

        end
    end
end

function plotSavedLines()
    global saveLines
    
    for i = 1:length(saveLines)
        nexttile(saveLines{i}{1});
        x = saveLines{i}{2};
        y = saveLines{i}{3};
        color = saveLines{i}{4};
        width = saveLines{i}{5};
        style = saveLines{i}{6};
        hold on;
        plot(x,y,'Color',color,"LineWidth",width,'LineStyle',style);
        hold off;
    end
end

function finalStyling()
    global axs minT maxT mode errFloor legendtile labelList
    nexttile(1);
    ylabel("Position Proxy, $\mathrm{q}_0 = \cos(\theta/2)$","FontSize",axs);
%     xlabel("Time (s)","FontSize",axs);
    xlim([minT,maxT]);
    ylim([-1  ,1]);
    grid on;
    set(gca,'xscale',mode);
    
    nexttile(2);
    ylabel("Absolute Angle Error,  $\left|\mathrm{q}_0 - \mathrm{q}_0^{true}\right|$","FontSize",axs);
    xlim([minT,maxT]);
    set(gca,'yscale','log');
%     xlabel("Time (s)","FontSize",axs);
    ylim(10.^[errFloor,0.5]);
    grid on;
    set(gca,'xscale',mode);
    
    nexttile(3);
    xlabel("Time (s)","FontSize",axs);
    ylabel("Relative Energy Error: $\left|\frac{E - E_{true}}{E_{true}}\right|$","FontSize",axs);
    xlim([minT,maxT]);
    set(gca,'yscale','log');
    ylim([10.^errFloor,1]);
    grid on;
    set(gca,'xscale',mode);

    nexttile(legendtile);
    legend(labelList,"FontSize",18,"Location","northwest");

    
end

function name = extractName(fileName)
    types = ["Brute","Mag","Sym"];
    renameTypes = ["Linear","Mag","Sym"];
    
    n = split(fileName,"_");
    i = contains(types,n(1));
    name = "\texttt{" + renameTypes(i) + "}";
    
    if contains(name,["Mag","Sym"])
        name = name + n(3) + " $\big($" + n(2) + ", $";
    else
        name = name + " $\big(";
    end
        

    res = str2num(extractBetween(n(end),2,strlength(n(end)) - 4));
    
    pow = floor(res/10);
    pref = res - 10*pow;
    s = "";
    if pref > 0
        s = s + num2str(round(10.^(pref/10))) + "\times ";
    end
    s = s + "10^{" + num2str(pow) + "}\big)$";

    name = name + s;
    

end