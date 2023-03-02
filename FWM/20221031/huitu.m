clc;
close all;
clear;

dirs = dir('*.xlsx');
N = length(dirs);

GC_220_1 = xlsread('1-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_2 = xlsread('1-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_3 = xlsread('1-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_4 = xlsread('1-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_5 = xlsread('1-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_6 = xlsread('1-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_7 = xlsread('1-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_8 = xlsread('1-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_9 = xlsread('1-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_10 = xlsread('1-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_11 = xlsread('1-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_16 = xlsread('2-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_17 = xlsread('2-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_18 = xlsread('2-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_19 = xlsread('2-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_20 = xlsread('2-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_21 = xlsread('2-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_22 = xlsread('2-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_23 = xlsread('2-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_24 = xlsread('2-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_25 = xlsread('2-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_26 = xlsread('2-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_31 = xlsread('3-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_32 = xlsread('3-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_33 = xlsread('3-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_34 = xlsread('3-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_35 = xlsread('3-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_36 = xlsread('3-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_37 = xlsread('3-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_38 = xlsread('3-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_39 = xlsread('3-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_40 = xlsread('3-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_41 = xlsread('3-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_45 = xlsread('4-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_46 = xlsread('4-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_47 = xlsread('4-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_48 = xlsread('4-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_49 = xlsread('4-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_50 = xlsread('4-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_51 = xlsread('4-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_52 = xlsread('4-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_53 = xlsread('4-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_54 = xlsread('4-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_55 = xlsread('4-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_56 = xlsread('5-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_57 = xlsread('5-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_58 = xlsread('5-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_59 = xlsread('5-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_60 = xlsread('5-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_61 = xlsread('5-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_62 = xlsread('5-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_63 = xlsread('5-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_64 = xlsread('5-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_65 = xlsread('5-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_66 = xlsread('5-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_68 = xlsread('6-1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_69 = xlsread('6-2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_70 = xlsread('6-3.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_71 = xlsread('6-4.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_72 = xlsread('6-5.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_73 = xlsread('6-6.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_74 = xlsread('6-7.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_75 = xlsread('6-8.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_76 = xlsread('6-9.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_77 = xlsread('6-10.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_78 = xlsread('6-11.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据

GC_220_79 = xlsread('1.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据
GC_220_80 = xlsread('2.xlsx', 'IL', 'A:B');  %读入220nm厚SOI的光栅耦合器数据




for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('1-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('1-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('1-%d.jpg',i));

%  end
end

for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('2-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('2-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('2-%d.jpg',i));

%  end
end

for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('3-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('3-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('3-%d.jpg',i));

%  end
end

for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('4-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('4-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('4-%d.jpg',i));

%  end
end

for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('5-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('5-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('5-%d.jpg',i));

%  end
end

for i=1:11
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('6-%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('6-%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('6-%d.jpg',i));

%  end
end

for i=1:2
%       for j= 1:9
          figure;  hold on;
file_name=sprintf('%g.xlsx',i);  %器件文件名称
%file_name=sprintf('cell2_TE1_bus_%g.xlsx',i);
path=strcat('D:\测试\AIN\21.3.31AIN Ring\涂胶\',file_name);
I=xls(path);
%I2=[I(:,1) I(:,2)+GC_220_80(:,2) ];
I2=[I(:,1) I(:,2)];
h = plot(I2(:,1),I2(:,2));

set(h, 'MarkerFaceColor', h.Color, 'linewidth', 2, 'markersize', 6);
axis([1507 1620 -70 -20])
xlabel('Wavelength (nm)')
ylabel('Transmission (dB)')
legend({sprintf('%d',i)}, 'location', 'southeast');
Plot_Setting();
saveas(gcf,sprintf('%d.jpg',i));

%  end
end
