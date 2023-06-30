%% Satellite Animation

function ChasePlot(xt,xc)

R = norm(xc(:,1))*2;

we = 2*pi/(24*3600);

load Data/EARTH_1;
load Data/EARTH_2;
load Data/EARTH_05;
load Data/COAST_3D;
load Data/COAST_C_3D;
load Data/COAST_L_3D;
load('Data/coast.mat');

theta = 0;
    
rECEFt = rotz(theta*180/pi)*xt;
rECEFc = rotz(theta*180/pi)*xc;

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
plot3(xt(1,:),xt(2,:),xt(3,:),'Color',[1 0 0],'Linewidth',2);
plot3(xc(1,:),xc(2,:),xc(3,:),'Color',[0 1 0],'Linewidth',2);
axis equal
axis off
grid off
box off
rotate3d on
axis([-R R -R R -R R]);
zoom(1.8);
hold off

subplot('Position',[0 0 1 .3]);%subplot(1,3,3); 
hold on;
%title('\color{white}Latitudinal-Longitudinal Plot');
pcolor(LON,LAT,TOPOGRAPHY);
demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
plot(longt,latt,'Color',[1 0 0],'Linewidth',2);
plot(longc,latc,'Color',[0 1 0],'Linewidth',2);
shading interp
axis equal
axis off
grid off
box off
hold off

end