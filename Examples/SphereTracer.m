set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',23)

colormap(gray)
global ell
ell = 0.06;


% files = "Output/Frog" + ["0"]+ "/Sym_Frog_2_N50.dat";
files = "Output/Frog0/" + ["Sym_Frog_2_N63.dat","Brute_N43.dat"];
styles = ["-","--",":"];
c1 = [0.2,1,0.2];
c2 = [0.8,0.2,0.2];
clf;
plotSphere();
tSet = {};
xSet = {};
zSet = {};

for i = 1:length(files)
    
    [x,z,t] = getFrameData(files(i));
    xSet{end+1} = x;
    zSet{end+1} = z;
    tSet{end+1} = t;
%     plotFrameN(x,z,round(length(x)*0.7),styles(i),c1*0.8^i,c2*0.8^i);
end
disp("Files read to memory");
mod = 0.8.^[1:length(xSet)];

c1s = [];
c2s = [];
for i = 1:length(xSet)
   c1s(end+1,:) = c1 * 0.8^i;
   c2s(end+1,:) = c2 * 0.8^i;
end

plotGroupFrames(xSet,zSet,tSet,1,1e7,30,styles,c1s,c2s,"TestVideo");

function [xp,zp,tp] = getFrameData(fileName)
    x = [0;1;0;0];
    z = [0;0;0;1];
    xp = [];
    zp = [];
    tp = [];
    f= readtable(fileName);

    for i = 1:1:height(f)
        tp(end+1) = f.t(i);
        q0 = f.q0(i);
        q1 = f.q1(i);
        q2 = f.q2(i);
        q3 = f.q3(i);



        M = Left(q0,q1,q2,q3) * Right(q0,-q1,-q2,-q3);

        Mx = M*x;
        Mz = M*z;
        xp(end+1,:) = Mx(2:end);
        zp(end+1,:) = Mz(2:end);
    end
    global ell
    xp = zp + ell*xp;
end
function plotSphere()
    [sx,sy,sz] = sphere(25);
    surf(sx,sy,sz,'LineStyle','None','FaceAlpha',0.2);
end
function plotFrameN(xp,zp,n,style,c1,c2)
    
    memory = 2e3;



    for i =  n:n
    
        zz = zp(i,:);
        xx = xp(i,:);
        lookback = max(1,i-memory);

    %     clf;
        hold on;
        width = 2;
        plot3([0,zz(1),xx(1)],[0,zz(2),xx(2)],[0,zz(3),xx(3)],'Color','b','LineWidth',2*width,'LineStyle',style);

        plot3(xp(lookback:i,1),xp(lookback:i,2),xp(lookback:i,3),'Color',[c2,0.3],'LineWidth',width,'LineStyle',style)
        plot3(zp(lookback:i,1),zp(lookback:i,2),zp(lookback:i,3),'Color',c1,'LineWidth',width,'LineStyle',style)

        hold off;
        
    end
    
end
function plotGroupFrames(xpSet,zpSet,tSet,n1,n2,gap,styles,c1s,c2s,fileName)
    
    video =VideoWriter(fileName);
    video.FrameRate = 10;
    open(video);
    [sx,sy,sz] = sphere(10  );
    
    for i = n1:gap:n2
        i
        clf;
%         surf(sx,sy,sz,'LineStyle','None','FaceAlpha',0.2);
        tit = "$ T = " + num2str(tSet{1}(i)) + "$ seconds";
        title(tit);
        for j = 1:length(xpSet)
            if (i < length(xpSet{j}))
               
                 plotFrameN(xpSet{j},zpSet{j},i,styles(j),c1s(j,:),c2s(j,:));
            end
        end
        view(45,27);
        xlim([-1,1]);
        ylim([-1,1]);
        zlim([-1,1]);
%         drawnow;
        writeVideo(video,getframe(gcf));
    end
    
% axis equal;
end
function L = Left(q0,q1,q2,q3)
    L = [ [q0,-q1,-q2,-q3];[q1,q0,-q3,q2];[q2,q3,q0,-q1];[q3,-q2,q1,q0]];
end
function R = Right(q0,q1,q2,q3)
     R = [ [q0,-q1,-q2,-q3];[q1,q0,q3,-q2];[q2,-q3,q0,q1];[q3,q2,-q1,q0]];
end