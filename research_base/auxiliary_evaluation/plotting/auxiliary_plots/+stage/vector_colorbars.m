function vector_colorbars(fig)

drawnow;bars=findall(fig,'Type','ColorBar');

for cb=reshape(bars,1,[])
 if strcmp(cb.TickLabelsMode,'auto')
  cb.TickLabels=cellstr(compose('%.3g',cb.Ticks));
  try,cb.Ruler.Exponent=0;catch,end
 end
end
drawnow;
for cb=reshape(bars,1,[])
 pos=cb.Position;lim=cb.Limits;ticks=cb.Ticks;labels=cb.TickLabels;label=cb.Label.String;
 if strcmp(cb.TickLabelsMode,'auto'),labels=cellstr(compose('%.3g',ticks));end
 try,target=cb.Axes;map=colormap(target);catch
  axs=findall(fig,'Type','axes');dd=inf(numel(axs),1);for k=1:numel(axs),ap=axs(k).Position;dd(k)=abs(ap(2)+ap(4)/2-pos(2)-pos(4)/2)+abs(ap(1)+ap(3)-pos(1));end;[~,kk]=min(dd);map=colormap(axs(kk));
 end
 cb.Visible='off';cb.Label.Visible='off';n=size(map,1);
 ax=axes(fig,'Position',pos,'Tag','VectorColorbar','FontName','Times New Roman','FontSize',8,'LineWidth',.7,...
  'XLim',[0 1],'YLim',lim,'XTick',[],'YTick',ticks,'YTickLabel',labels,'YAxisLocation','right','Box','on','TickDir','out','Layer','top');hold(ax,'on');
 for j=1:n,y1=lim(1)+(j-1)/n*diff(lim);y2=lim(1)+j/n*diff(lim);patch(ax,[0 1 1 0],[y1 y1 y2 y2],map(j,:),'EdgeColor','none');end
 ylabel(ax,label,'FontName','Times New Roman','FontSize',9);
end
end
