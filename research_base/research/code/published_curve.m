function [x,y,hx,hy]=published_curve(T,caseid,series,regime,fraction)

if nargin<5,fraction=0;end
if string(regime)=="m",region=string(T.quantity)=="moment";else,region=endsWith(string(T.regime),regime);end
q=T(string(T.case_id)==string(caseid)&string(T.series)==string(series)&...
 region&T.load_fraction==fraction&T.visible==1,:);
assert(~isempty(q),'No matching published curve');
if string(q.quantity(1))=="moment",x=q.y_value;y=q.x_value;hx=q.y_halfstroke;hy=q.x_halfstroke;
else,x=q.x_value;y=q.y_value;hx=q.x_halfstroke;hy=q.y_halfstroke;end
ok=x>=0&y>=0;x=x(ok);y=y(ok);hx=hx(ok);hy=hy(ok);
[x,ia]=unique(x);y=y(ia);hx=hx(ia);hy=hy(ia);
if fraction==0 && x(1)>0,x=[0;x];y=[0;y];hx=[hx(1);hx];hy=[hy(1);hy];end
end
