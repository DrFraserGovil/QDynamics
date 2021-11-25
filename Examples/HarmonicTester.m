set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
figure(2);
axs = 22;

LEGENDPLOT = 1;

folder = "Output/Harmonic";
q = dir(folder);
clf;
T = tiledlayout(2,2,"TileSpacing","Compact","Padding","None");

minT = 0.1 ;
maxT = 30  ;
Nt = 2000;
t = 10.^linspace(log10(minT),log10(maxT),Nt);
initialAngle = 1;
U0 = 1;
Jx = 1;

angle = @(t) initialAngle * cos(sqrt(2*U0/(Jx))*t);
momentum = @(t) -initialAngle * sqrt(2*U0*Jx) * sin(sqrt(2*U0/Jx) * t);
potential = @(t) U0 * angle(t).^2;
energy = @(t) momentum(t).^2/(2* Jx) + potential(t);
analyticTheta = angle(t);

nexttile(1);
plot(t,cos(analyticTheta/2),'k');
nexttile(LEGENDPLOT);
plot(NaN,NaN,'k');
hold on;
% plot(t,(analyticTheta),'k',"HandleVisibility","off");
% plot(t,cummin(analyticTheta),'k',"HandleVisibility","off");
hold off;

fs = ["Analytic Solution"];
cols = [colororder; rand(14,3)];
line = 0;
width = 2   ;
global errFloor;
errFloor = -6;


triggerList = ["Brute","Mag_Euler_0","Mag_Euler_1","Mag_Euler_2","Sym_Euler_2","Sym_Frog"];
% triggerList = ["Brute","SymiL2"];
triggerId = zeros(size(triggerList))-1;
typeList = ["-","--",":"];
for f = {q.name}
	file = f{1};
    
    trig = -1;
    for i = 1:length(triggerList)
        if contains(file,triggerList(i))
            trig = i;
            triggerId(i) = triggerId(i) + 1;
        end
    end
    
	if trig > 0
		color = cols(trig,:);
        sID = mod(triggerId(trig),length(typeList))+1;
        loopies = floor(triggerId(trig)/length(typeList));
        color = color * 0.7^loopies;
        style = typeList(sID);
        nexttile(LEGENDPLOT)
%         plot(NaN,NaN,'Color',color,'LineStyle',style,'LineWidth',width);
        line = line + 1;
        
		f = readtable(folder + "/" +file);	
        if height(f) > 0
            file

            fs(end+1) = "\verb|" + file+ "|";
            theta = f.q0;
            theta(isnan(theta)) = 1e-15;

            nexttile(1);
            hold on;
    %         plot(f.t,(theta),'Color',cols(line,:),"LineWidth",width);

            envelopePlot(f.t,theta,color,width,style);
            hold off;
    % 		xlabel("Time (s)","FontSize",axs);
            ylabel("Rotation Angle, $\cos(\theta/2)$","FontSize",axs);

            nexttile(2);
            hold on;
            orig = angle(f.t);
            diff = abs(theta - orig);
            smoothPlot(f.t,diff,color,width,style);
            hold off;
    % 		xlabel("Time (s)","FontSize",axs);
            ylabel("Angle Error, $\log_{10}\left(|\theta - \theta_{true}|\right)$","FontSize",axs);

            nexttile(3);
            hold on;
            orig = 0;
            diff = abs(f.H - orig);
            smoothPlot(f.t,diff,color,width,style);
            hold off;
            xlabel("Time (s)","FontSize",axs);
            ylabel("Energy Error, $\log_{10}\left(|E - E{true}|\right)$","FontSize",axs);

            nexttile(4);
            hold on;
            orig = 0;%momentum(f.t);
            diff = abs(abs(f.L) - abs(orig));

            smoothPlot(f.t,diff,color,width,style);
            hold off;
            xlabel("Time (s)","FontSize",axs);
            ylabel("Momentum Error, $\log_{10}\left(|L - L{true}|\right)$","FontSize",axs);
%             pause(3)
            drawnow;
        end
        
    end
end
mode = 'log';
nexttile(1);
xlim([min(t),max(t)]);
ylim([-1   ,1]);
set(gca,'xscale',mode);
nexttile(2);



xlim([min(t),max(t)]);
ylim([errFloor,0.6]);
set(gca,'xscale',mode);
nexttile(3);
xlim([min(t),max(t)]);
% ylim([0,20   ]);
set(gca,'xscale',mode);
nexttile(4);
xlim([min(t),max(t)]);
ylim([0,3  ]);
set(gca,'xscale',mode);

nexttile(LEGENDPLOT)
legend(fs,"FontSize",12,"Location","southeast");


nexttile(1);
hold on;
plot(t,cos(analyticTheta/2),'k');
hold off;
nexttile(3);
hold on;
plot(t,energy(t))
plot(t,potential(t));
hold off;

nexttile(4);
hold on;
plot(t,abs(momentum(t)));
hold off;
function y =clarifyer(x)
    global errFloor;
    errorMag = 20000;
    
    err = 10^(errFloor-1);

    y = (x + err)*(1.0+rand/errorMag);
end
function y = smoother(x)
    y = smooth(x,10);
%     y = cummax(x);
end
function envelopePlot(x,y,color,width,style)

%     sampler =100*ones(size(y));
%     [~,wp] = findpeaks(y);
%     [~,wt] = findpeaks(-y);
%     xp = x(wp);
%     yp = y(wp);
%     xt = x(wt);
%     yt = y(wt);

    plot(x,y,'Color',color,"LineWidth",width,'LineStyle',style);
%     plot(xp,yp,'Color',color,"LineWidth",width,'LineStyle',style);
% 	plot(xt,yt,'Color',color,"LineWidth",width,"HandleVisibility","off",'LineStyle',style);
end
function smoothPlot(x,diff,color,width,style)
    y = (clarifyer(diff));

    smoothing = 50;
    goodPoints = sum(~isnan(y));
    
    while goodPoints < smoothing
        smoothing = smoothing/5;
    end

    smoothed = smooth(y,smoothing);
    
    plot(x(1:end-smoothing),smoothed(1:end-smoothing),'Color',color,"LineWidth",width,'LineStyle',style);
   
end