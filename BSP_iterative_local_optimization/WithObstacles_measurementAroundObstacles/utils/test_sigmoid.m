clc;
clear all
close all

x = [0:0.1:100];
for i = 1:length(x)
    y(i) = 1/(1+exp((x(i)-50))) - 0.5;
    z(i) = 1/( 1+exp( -1*(x(i)-70) ) ) - 0.5;
   
end
figure;
plot(x,y);
hold on;
plot(x,z);

% figure;
% mapx = 0:100;
% mapy = 0:100;
% [mapX,mapY] = meshgrid(mapx,mapy);
% 
% % contZvar = Amp*exp(-(mapX - mu).^2/sig);
% 
% for i = 1:length(mapx)
%     for j = 1 : length(mapy)
%         contZvar(i,j) = (1/(1+exp(-(mapx(i)-50))) + 1/(1+exp((mapx(i)-70))) - 1)*(1/(1+exp(-(mapy(j)-50))) + 1/(1+exp((mapy(j)-70))) - 1);
%        
%     end
% end
% conto = contourf(mapX,mapY,contZvar,50, 'edgecolor','none');
% colormap((gray));