set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',23)

colormap(gray)
global ell
ell = 0.06;

files = "Output/Frog" + ["0","1","2"]+ "/Sym_Frog_2_N50.dat";
styles = ["-","--",":"];
c1 = [0.2,1,0.2];
c2 = [0.8,0.2,0.2];
clf;
plotSphere();
for i = 1:length(files)
    
    [x,z] = getFrameData(files(i));
    plotFrameN(x,z,round(length(x)*0.7),styles(i),c1*0.8^i,c2*0.8^i);
end
view(45,7);
axis equal;
function [xp,zp] = getFrameData(fileName)
    x = [0;1;0;0];
    z = [0;0;0;1];
    xp = [];
    zp = [];
    f= readtable(fileName);

    for i = 1:1:height(f)

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
    
    memory = 1e5;



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

function L = Left(q0,q1,q2,q3)
    L = [ [q0,-q1,-q2,-q3];[q1,q0,-q3,q2];[q2,q3,q0,-q1];[q3,-q2,q1,q0]];
end
function R = Right(q0,q1,q2,q3)
     R = [ [q0,-q1,-q2,-q3];[q1,q0,q3,-q2];[q2,-q3,q0,q1];[q3,q2,-q1,q0]];
end