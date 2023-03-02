close all
clear all
hold on


for i = 145
    if i ~= 125
        if i <100
            str = sprintf('%s%s%s','W00',num2str(i),'.CSV');
        else
            str = sprintf('%s%s%s','W0',num2str(i),'.CSV');
        end
    end


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
 xlim([1400, 1650])     
ylim([-100, 0])

hhh = string;
for i = 1:(7)
    hhh(i) =sprintf('%s%s','Signal wavelength:',num2str(i*2+1550),'nm') ;
end
% legend(hhh);
% 
