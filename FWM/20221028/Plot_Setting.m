function Plot_Setting()
    grid on; box on;
    set(gca, 'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 20, 'GridLineStyle', '--');
    width=1000;%��ȣ�������
    height=500;%�߶�
    left=200;%����Ļ���½�ˮƽ����
    bottem=100;%����Ļ���½Ǵ�ֱ����
    set(gcf,'position',[left,bottem,width,height]);
    
end