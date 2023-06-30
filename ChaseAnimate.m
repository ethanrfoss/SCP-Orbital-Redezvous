%% Satellite Animation

function ChaseAnimate(tr,xt,xc,P,saveName)

dt = P/100; % Sample rate of data [s]

t = 0:dt:tr(end);

rt = interp1(tr,xt',t)';
rc = interp1(tr,xc',t)';

R = norm(rc(:,1))*2;

we = 2*pi/(24*3600);

trail = .5*ceil(P/dt);
traili = ceil(linspace(1,trail,6));

load Data/EARTH_1;
load Data/EARTH_2;
load Data/EARTH_05;
load Data/COAST_3D;
load Data/COAST_C_3D;
load Data/COAST_L_3D;
load('Data/coast.mat');

theta = t*we;

for i = 1:length(t)
    
    rECEFt(:,i) = rotz(theta(i)*180/pi)*rt(:,i);
    rECEFc(:,i) = rotz(theta(i)*180/pi)*rc(:,i);
    
end

llat = ecef2lla(rECEFt')';
latt = llat(1,:); longt = mod(llat(2,:)+360,360);
longt(abs(diff(longt))>90) = NaN;
llac = ecef2lla(rECEFc')';
latc = llac(1,:); longc = mod(llac(2,:)+360,360);
longc(abs(diff(longc))>90) = NaN;

f = figure('Position',[0 50 800 950]);
set(f,'Color','k');
subplot('Position',[0 .333 1 .6666]);%subplot(1,3,1); 
hold on; view(45,45);
%title('\color{white}Earth Rotating Frame');
pEarth = surf(X_EARTH,Y_EARTH,Z_EARTH,TOPOGRAPHY);
demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
pCoast = plot3(X_COAST,Y_COAST,Z_COAST,'Color','k');
shading interp;
% pECI = plot3(r(1,1:2),r(2,1:2),r(3,1:2),'r','Linewidth',2);
pECI1t = plot3(rt(1,traili(1):traili(2)),rt(2,traili(1):traili(2)),rt(3,traili(1):traili(2)),'Color',[1 0 0 .2],'Linewidth',2);
pECI2t = plot3(rt(1,traili(2):traili(3)),rt(2,traili(2):traili(3)),rt(3,traili(2):traili(3)),'Color',[1 0 0 .4],'Linewidth',2);
pECI3t = plot3(rt(1,traili(3):traili(4)),rt(2,traili(3):traili(4)),rt(3,traili(3):traili(4)),'Color',[1 0 0 .6],'Linewidth',2);
pECI4t = plot3(rt(1,traili(4):traili(5)),rt(2,traili(4):traili(5)),rt(3,traili(4):traili(5)),'Color',[1 0 0 .8],'Linewidth',2);
pECI5t = plot3(rt(1,traili(5):traili(6)),rt(2,traili(5):traili(6)),rt(3,traili(5):traili(6)),'Color',[1 0 0 1],'Linewidth',2);
pECI1c = plot3(rc(1,traili(1):traili(2)),rc(2,traili(1):traili(2)),rc(3,traili(1):traili(2)),'Color',[0 1 0 .2],'Linewidth',2);
pECI2c = plot3(rc(1,traili(2):traili(3)),rc(2,traili(2):traili(3)),rc(3,traili(2):traili(3)),'Color',[0 1 0 .4],'Linewidth',2);
pECI3c = plot3(rc(1,traili(3):traili(4)),rc(2,traili(3):traili(4)),rc(3,traili(3):traili(4)),'Color',[0 1 0 .6],'Linewidth',2);
pECI4c = plot3(rc(1,traili(4):traili(5)),rc(2,traili(4):traili(5)),rc(3,traili(4):traili(5)),'Color',[0 1 0 .8],'Linewidth',2);
pECI5c = plot3(rc(1,traili(5):traili(6)),rc(2,traili(5):traili(6)),rc(3,traili(5):traili(6)),'Color',[0 1 0 1],'Linewidth',2);
sECIt = scatter3(rt(1,traili(6)),rt(2,traili(6)),rt(3,traili(6)),60,'w','filled');
sECIc = scatter3(rc(1,traili(6)),rc(2,traili(6)),rc(3,traili(6)),60,'w','filled');
axis equal
axis off
grid off
box off
rotate3d on
axis([-R R -R R -R R]);
zoom(1.8);
hold off
% subplot('Position',[.5 .333 .5 .6666]);%subplot(1,3,2); 
% hold on; view(45,45);
% %title('\color{white}Earth Fixed Frame');
% surf(X_EARTH,Y_EARTH,Z_EARTH,TOPOGRAPHY);
% demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
% plot3(X_COAST,Y_COAST,Z_COAST,'Color','k');
% %pECEF = plot3(rECEF(1,1),rECEF(2,1),rECEF(3,1),'r','Linewidth',2);
% pECEF1t = plot3(rECEFt(1,traili(1):traili(2)),rECEFt(2,traili(1):traili(2)),rECEFt(3,traili(1):traili(2)),'Color',[1 0 0 .2],'Linewidth',2);
% pECEF2t = plot3(rECEFt(1,traili(2):traili(3)),rECEFt(2,traili(2):traili(3)),rECEFt(3,traili(2):traili(3)),'Color',[1 0 0 .4],'Linewidth',2);
% pECEF3t = plot3(rECEFt(1,traili(3):traili(4)),rECEFt(2,traili(3):traili(4)),rECEFt(3,traili(3):traili(4)),'Color',[1 0 0 .6],'Linewidth',2);
% pECEF4t = plot3(rECEFt(1,traili(4):traili(5)),rECEFt(2,traili(4):traili(5)),rECEFt(3,traili(4):traili(5)),'Color',[1 0 0 .8],'Linewidth',2);
% pECEF5t = plot3(rECEFt(1,traili(5):traili(6)),rECEFt(2,traili(5):traili(6)),rECEFt(3,traili(5):traili(6)),'Color',[1 0 0 1],'Linewidth',2);
% pECEF1c = plot3(rECEFc(1,traili(1):traili(2)),rECEFc(2,traili(1):traili(2)),rECEFc(3,traili(1):traili(2)),'Color',[0 1 0 .2],'Linewidth',2);
% pECEF2c = plot3(rECEFc(1,traili(2):traili(3)),rECEFc(2,traili(2):traili(3)),rECEFc(3,traili(2):traili(3)),'Color',[0 1 0 .4],'Linewidth',2);
% pECEF3c = plot3(rECEFc(1,traili(3):traili(4)),rECEFc(2,traili(3):traili(4)),rECEFc(3,traili(3):traili(4)),'Color',[0 1 0 .6],'Linewidth',2);
% pECEF4c = plot3(rECEFc(1,traili(4):traili(5)),rECEFc(2,traili(4):traili(5)),rECEFc(3,traili(4):traili(5)),'Color',[0 1 0 .8],'Linewidth',2);
% pECEF5c = plot3(rECEFc(1,traili(5):traili(6)),rECEFc(2,traili(5):traili(6)),rECEFc(3,traili(5):traili(6)),'Color',[0 1 0 1],'Linewidth',2);
% sECEFt = scatter3(rECEFt(1,traili(6)),rECEFt(2,traili(6)),rECEFt(3,traili(6)),60,'w','filled');
% sECEFc = scatter3(rECEFc(1,traili(6)),rECEFc(2,traili(6)),rECEFc(3,traili(6)),60,'w','filled');
% shading interp;
% axis equal
% axis off
% grid off
% box off
% rotate3d on
% axis([-R R -R R -R R]);
% zoom(1.8);
% hold off
subplot('Position',[0 0 1 .3]);%subplot(1,3,3); 
hold on;
%title('\color{white}Latitudinal-Longitudinal Plot');
pcolor(LON,LAT,TOPOGRAPHY);
demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
%pLL = plot(long(1:2),lat(1:2),'r','Linewidth',2);
pLL1t = plot(longt(traili(1):traili(2)),latt(traili(1):traili(2)),'Color',[1 0 0 .2],'Linewidth',2);
pLL2t = plot(longt(traili(2):traili(3)),latt(traili(2):traili(3)),'Color',[1 0 0 .4],'Linewidth',2);
pLL3t = plot(longt(traili(3):traili(4)),latt(traili(3):traili(4)),'Color',[1 0 0 .6],'Linewidth',2);
pLL4t = plot(longt(traili(4):traili(5)),latt(traili(4):traili(5)),'Color',[1 0 0 .8],'Linewidth',2);
pLL5t = plot(longt(traili(5):traili(6)),latc(traili(5):traili(6)),'Color',[1 0 0 1],'Linewidth',2);
pLL1c = plot(longc(traili(1):traili(2)),latc(traili(1):traili(2)),'Color',[0 1 0 .2],'Linewidth',2);
pLL2c = plot(longc(traili(2):traili(3)),latc(traili(2):traili(3)),'Color',[0 1 0 .4],'Linewidth',2);
pLL3c = plot(longc(traili(3):traili(4)),latc(traili(3):traili(4)),'Color',[0 1 0 .6],'Linewidth',2);
pLL4c = plot(longc(traili(4):traili(5)),latc(traili(4):traili(5)),'Color',[0 1 0 .8],'Linewidth',2);
pLL5c = plot(longc(traili(5):traili(6)),latc(traili(5):traili(6)),'Color',[0 1 0 1],'Linewidth',2);
sLLt = scatter(longt(traili(6)),latt(traili(6)),60,'w','filled');
sLLc = scatter(longc(traili(6)),latc(traili(6)),60,'w','filled');
shading interp
axis equal
axis off
grid off
box off
hold off

[row col] = size(X_EARTH);

POS_EARTH = [X_EARTH(:)'; Y_EARTH(:)'; Z_EARTH(:)'];
POS_COAST = [X_COAST';Y_COAST';Z_COAST'];

step = 1;
for i = 1+step:step:length(t)-traili(6)
    
    if ishandle(f) == 0
        break;
    end
    
    Earth = rotz(-theta(i)*180/pi)*POS_EARTH;
    X_EARTH = reshape(Earth(1,:),row,col);
    Y_EARTH = reshape(Earth(2,:),row,col);
    Z_EARTH = reshape(Earth(3,:),row,col);
    Coast = rotz(-theta(i)*180/pi)*POS_COAST;
    
    set(pEarth,'XData',X_EARTH,'YData',Y_EARTH,'ZData',Z_EARTH);
    set(pCoast,'XData',Coast(1,:),'YData',Coast(2,:),'ZData',Coast(3,:));
    
    set(sECIt,'XData',rt(1,traili(6)+i),'YData',rt(2,traili(6)+i),'ZData',rt(3,traili(6)+i));
    set(sECIc,'XData',rc(1,traili(6)+i),'YData',rc(2,traili(6)+i),'ZData',rc(3,traili(6)+i));
%     set(sECEFt,'XData',rECEFt(1,traili(6)+i),'YData',rECEFt(2,traili(6)+i),'ZData',rECEFt(3,traili(6)+i));
%     set(sECEFc,'XData',rECEFc(1,traili(6)+i),'YData',rECEFc(2,traili(6)+i),'ZData',rECEFc(3,traili(6)+i));
    set(sLLt,'XData',longt(traili(6)+i),'YData',latt(traili(6)+i));
    set(sLLc,'XData',longc(traili(6)+i),'YData',latc(traili(6)+i));
    
    %set(pECI,'XData',r(1,max(i-trail,1):i),'YData',r(2,max(i-trail,1):i),'ZData',r(3,max(i-trail,1):i));
    set(pECI1t,'XData',rt(1,[traili(1):traili(2)]+i),'YData',rt(2,[traili(1):traili(2)]+i),'ZData',rt(3,[traili(1):traili(2)]+i));
    set(pECI2t,'XData',rt(1,[traili(2):traili(3)]+i),'YData',rt(2,[traili(2):traili(3)]+i),'ZData',rt(3,[traili(2):traili(3)]+i));
    set(pECI3t,'XData',rt(1,[traili(3):traili(4)]+i),'YData',rt(2,[traili(3):traili(4)]+i),'ZData',rt(3,[traili(3):traili(4)]+i));
    set(pECI4t,'XData',rt(1,[traili(4):traili(5)]+i),'YData',rt(2,[traili(4):traili(5)]+i),'ZData',rt(3,[traili(4):traili(5)]+i));
    set(pECI5t,'XData',rt(1,[traili(5):traili(6)]+i),'YData',rt(2,[traili(5):traili(6)]+i),'ZData',rt(3,[traili(5):traili(6)]+i));
    set(pECI1c,'XData',rc(1,[traili(1):traili(2)]+i),'YData',rc(2,[traili(1):traili(2)]+i),'ZData',rc(3,[traili(1):traili(2)]+i));
    set(pECI2c,'XData',rc(1,[traili(2):traili(3)]+i),'YData',rc(2,[traili(2):traili(3)]+i),'ZData',rc(3,[traili(2):traili(3)]+i));
    set(pECI3c,'XData',rc(1,[traili(3):traili(4)]+i),'YData',rc(2,[traili(3):traili(4)]+i),'ZData',rc(3,[traili(3):traili(4)]+i));
    set(pECI4c,'XData',rc(1,[traili(4):traili(5)]+i),'YData',rc(2,[traili(4):traili(5)]+i),'ZData',rc(3,[traili(4):traili(5)]+i));
    set(pECI5c,'XData',rc(1,[traili(5):traili(6)]+i),'YData',rc(2,[traili(5):traili(6)]+i),'ZData',rc(3,[traili(5):traili(6)]+i));
    %set(pECEF,'XData',rECEF(1,max(i-trail,1):i),'YData',rECEF(2,max(i-trail,1):i),'ZData',rECEF(3,max(i-trail,1):i));
%     set(pECEF1t,'XData',rECEFt(1,[traili(1):traili(2)]+i),'YData',rECEFt(2,[traili(1):traili(2)]+i),'ZData',rECEFt(3,[traili(1):traili(2)]+i));
%     set(pECEF2t,'XData',rECEFt(1,[traili(2):traili(3)]+i),'YData',rECEFt(2,[traili(2):traili(3)]+i),'ZData',rECEFt(3,[traili(2):traili(3)]+i));
%     set(pECEF3t,'XData',rECEFt(1,[traili(3):traili(4)]+i),'YData',rECEFt(2,[traili(3):traili(4)]+i),'ZData',rECEFt(3,[traili(3):traili(4)]+i));
%     set(pECEF4t,'XData',rECEFt(1,[traili(4):traili(5)]+i),'YData',rECEFt(2,[traili(4):traili(5)]+i),'ZData',rECEFt(3,[traili(4):traili(5)]+i));
%     set(pECEF5t,'XData',rECEFt(1,[traili(5):traili(6)]+i),'YData',rECEFt(2,[traili(5):traili(6)]+i),'ZData',rECEFt(3,[traili(5):traili(6)]+i));
%     set(pECEF1c,'XData',rECEFc(1,[traili(1):traili(2)]+i),'YData',rECEFc(2,[traili(1):traili(2)]+i),'ZData',rECEFc(3,[traili(1):traili(2)]+i));
%     set(pECEF2c,'XData',rECEFc(1,[traili(2):traili(3)]+i),'YData',rECEFc(2,[traili(2):traili(3)]+i),'ZData',rECEFc(3,[traili(2):traili(3)]+i));
%     set(pECEF3c,'XData',rECEFc(1,[traili(3):traili(4)]+i),'YData',rECEFc(2,[traili(3):traili(4)]+i),'ZData',rECEFc(3,[traili(3):traili(4)]+i));
%     set(pECEF4c,'XData',rECEFc(1,[traili(4):traili(5)]+i),'YData',rECEFc(2,[traili(4):traili(5)]+i),'ZData',rECEFc(3,[traili(4):traili(5)]+i));
%     set(pECEF5c,'XData',rECEFc(1,[traili(5):traili(6)]+i),'YData',rECEFc(2,[traili(5):traili(6)]+i),'ZData',rECEFc(3,[traili(5):traili(6)]+i));
    %set(pLL,'XData',long(max(i-trail,1):i),'YData',lat(max(i-trail,1):i));
    set(pLL1t,'XData',longt([traili(1):traili(2)]+i),'YData',latt([traili(1):traili(2)]+i));
    set(pLL2t,'XData',longt([traili(2):traili(3)]+i),'YData',latt([traili(2):traili(3)]+i));
    set(pLL3t,'XData',longt([traili(3):traili(4)]+i),'YData',latt([traili(3):traili(4)]+i));
    set(pLL4t,'XData',longt([traili(4):traili(5)]+i),'YData',latt([traili(4):traili(5)]+i));
    set(pLL5t,'XData',longt([traili(5):traili(6)]+i),'YData',latt([traili(5):traili(6)]+i));
    set(pLL1c,'XData',longc([traili(1):traili(2)]+i),'YData',latc([traili(1):traili(2)]+i));
    set(pLL2c,'XData',longc([traili(2):traili(3)]+i),'YData',latc([traili(2):traili(3)]+i));
    set(pLL3c,'XData',longc([traili(3):traili(4)]+i),'YData',latc([traili(3):traili(4)]+i));
    set(pLL4c,'XData',longc([traili(4):traili(5)]+i),'YData',latc([traili(4):traili(5)]+i));
    set(pLL5c,'XData',longc([traili(5):traili(6)]+i),'YData',latc([traili(5):traili(6)]+i));
    
    drawnow;
    
    if nargin>4
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1+step
            imwrite(imind,cm,[cd '\' saveName '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
        else
            imwrite(imind,cm,[cd '\' saveName '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
        end
    end
    
end

end