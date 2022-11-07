function XTickRotateJMP(labs,rot)
    
    if nargin==1 || ~isnumeric(rot) || isinf(rot) || isnan(rot)
        rot=90;
    end
    
    n=length(labs);
    xl=xlim;
    fn=get(gca,'FontName');
    fs=get(gca,'FontSize');
    fw=get(gca,'FontWeight');
    xt=get(gca,'XTick');
    set(gca,'XTick',xt,'XTickLabel',[])
    if strcmp(get(gca,'YDir'),'reverse')
        for k=1:n
            text(xt(k),max(ylim),[labs{k},' '],'Rotation',rot,'HorizontalAlignment','right','VerticalAlignment','middle','FontName',fn,'FontSize',fs,'FontWeight',fw)
        end
    else
        for k=1:n
            text(xt(k),min(ylim),[labs{k},' '],'Rotation',rot,'HorizontalAlignment','right','VerticalAlignment','middle','FontName',fn,'FontSize',fs,'FontWeight',fw)
        end
    end
    xlim(xl)
    
end