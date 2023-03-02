close all
clear all
hold on


for i = 88:108
str = sprintf('%s%s%s','W00',num2str(i),'.CSV');    
data = csvread(str,30,0);
wavelength =data(:,1); 
loss1 =data(:,2); 
plot(wavelength,loss1,'linewidth',2);
end
% data = csvread('W0052.CSV'+i,30,0);
% wavelength =data(:,1); 
% loss2 =data(:,2); 
% plot(wavelength,loss2,'linewidth',2);
xlabel('Wavelength (nm)')
ylabel('Power (dBm)')
% legend("TM","TE");
Plot_Setting();
 xlim([1450, 1650])     
ylim([-100, 0])
hhh = string;
for i = 1:(87-74+1)
    hhh(i) =sprintf('%s%s','Signal wavelength:',num2str(i*5+1475),'nm') ;
end
legend(hhh);

