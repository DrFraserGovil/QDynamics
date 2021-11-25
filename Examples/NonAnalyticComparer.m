set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%% options
global errFloor legendtile folder width axs mode minT maxT;
errFloor = -10;
axs = 22;
legendtile = 1;
width = 2;
minT = 1e-1 ;
maxT = 30;
folder = "Output/Harmonic";
baseLineFile = "Sym_Frog_2_N51.dat";
triggerList = ["Sym_Euler_1","Sym_Euler_2"];
mode = 'log';
%% some persistent global quantities
global labelList triggerId cols positionInterp energyInterp momentumInterp
labelList = ["Baseline Model"];
triggerId = zeros(size(triggerList))-1;
cols = [colororder; rand(14,3)];
%% file checks
fileList = organiseFiles(folder,baseLineFile);

%% main plotting loop
preparePlot();
plotBase(baseLineFile)

for file = fileList
    plotFile(file,triggerList)
end

finalStyling();

function y =clarifyer(x)
    global errFloor;
    errorMag = 20000;
    
    err = 10^(errFloor-1);

    y = (x + err)*(1.0+rand/errorMag);
end

function envelopePlot(x,y,color,width,style,active)
    global errFloor minT
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
           s = smooth(y,length(y)/1000);
          plot(x,cummax(smooth(y,10)),'Color',color,"LineWidth",width,'LineStyle',style);
         
          plot(x,smooth(y,s),'Color',[color,0.1],"LineWidth",width,"HandleVisibility","off",'LineStyle',style);
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
    
    if sum(r==target) == 1
       s = target;
    else
        error("Target file, " + target + " could not be found, or there were multiple copies");
    end
    r(r == s) = [];
end
function preparePlot()  
    figure(1)
    clf;
    tiledlayout(2,2,"TileSpacing","Compact","Padding","Compact");
end
function plotBase(file)
     global legendtile folder positionInterp energyInterp momentumInterp errFloor;
     nexttile(1,[1,2]);
    baseW = 4;
    f = readtable(folder + "/" +file);	
    hold on;
    envelopePlot(f.t,f.q0,'k',baseW,'-',false);
    hold off;
    
    
    nexttile(4);
    hold on;
    yy = log10(abs(f.H - f.H(1)) + 10^(errFloor - 1));
    envelopePlot(f.t,yy,[0.1,0.1,0.1],baseW,'-',true);
    hold off;
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
    global triggerId folder cols labelList width positionInterp energyInterp momentumInterp
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
        if height(f) > 0
     
            labelList(end+1) = "\verb|" + fileName+ "|";
            theta = f.q0;
            theta(isnan(theta)) = 1e-15;

            nexttile(1,[1,2]);
            hold on;
    %         plot(f.t,(theta),'Color',cols(line,:),"LineWidth",width);

            envelopePlot(f.t,theta,color,width,style,false);
            hold off;
    % 		xlabel("Time (s)","FontSize",axs);
            

            nexttile(3);
            hold on;
            orig = positionInterp(f.t);
            diff = log10(abs(theta - orig));
            envelopePlot(f.t,diff,color,width,style,true);
            hold off;
    % 		xlabel("Time (s)","FontSize",axs);
            

            nexttile(4);
            hold on;
            orig = f.H(1);%energyInterp(f.t);
            diff = log10(abs(f.H - orig));
            envelopePlot(f.t,diff,color,width,style,true);
            hold off;
%            
% 
%             nexttile(4);
%             hold on;
%             orig = momentumInterp(f.t);
%             diff = log10(abs(abs(f.L) - abs(orig)));
% 
%             envelopePlot(f.t,diff,color,width,style,true);
%             hold off;
            
        end
    end
end

function finalStyling()
    global axs minT maxT mode errFloor legendtile labelList
    nexttile(1,[1,2]);
    ylabel("Position Proxy, $\mathrm{q}_0 = \cos(\theta/2)$","FontSize",axs);
    xlabel("Time (s)","FontSize",axs);
    xlim([minT,maxT]);
    ylim([-1   ,1]);
    grid on;
    set(gca,'xscale',mode);
    
    nexttile(3);
    ylabel("Angle Error, $\log_{10}\left(|\mathrm{q}_0 - \mathrm{q}_{0,true}|\right)$","FontSize",axs);
    xlim([minT,maxT]);
    xlabel("Time (s)","FontSize",axs);
%     ylim([errFloor,1.2]);
    grid on;
    set(gca,'xscale',mode);
    
    nexttile(4);
    xlabel("Time (s)","FontSize",axs);
    ylabel("Energy Error, $\log_{10}\left(|E - E{true}|\right)$","FontSize",axs);
    xlim([minT,maxT]);
% %     ylim([errFloor,1]);
    grid on;
    set(gca,'xscale',mode);
%     
%     nexttile(4);
%     xlabel("Time (s)","FontSize",axs);
%     ylabel("Momentum Error, $\log_{10}\left(|L - L{true}|\right)$","FontSize",axs);
%     xlim([minT,maxT]);
% %     ylim([errFloor,1]);
%     set(gca,'xscale',mode);
 
    nexttile(legendtile);
    legend(labelList,"FontSize",12,"Location","southwest");

    
end