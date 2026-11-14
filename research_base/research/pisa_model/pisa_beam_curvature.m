function curvature=pisa_beam_curvature(a,q)

ne=numel(a.z_m)-1;h=diff(a.z_m);v=q(1:2:2*(ne+1));t=q(2:2:2*(ne+1));g=q(2*(ne+1)+(1:ne));
eta=(v(1:end-1)-v(2:end))./h+(t(1:end-1)+t(2:end))/2+g;
dt=(t(2:end)-t(1:end-1))./h;
gp=[-.8611363115940526;-.3399810435848563;.3399810435848563;.8611363115940526];
curvature=reshape(6*gp*(eta./h)'+ones(4,1)*dt',[],1);
end
