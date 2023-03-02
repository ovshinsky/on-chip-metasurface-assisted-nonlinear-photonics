function Plot_Setting()
    grid on; box on;
    set(gca, 'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 20, 'GridLineStyle', '--');
    width=1000;%宽度，像素数
    height=500;%高度
    left=200;%距屏幕左下角水平距离
    bottem=100;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height]);
    
end