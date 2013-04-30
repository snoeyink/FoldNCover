function [result,fourpts] = enclsquaresize(p)
% p is hull vertices

    function result = perp(u)
        result = [-u(2) u(1)];
    end

result = Inf;
idx=convhull(p);
n = size(idx);
p = p([idx;idx(2:end)],:);
ix = [2 0 0 0];
u=p(2,:)-p(1,:);% initial direction
ln=norm(u);
v = [-perp(u);u;perp(u);-u]'; % directions for extremes
for j = 2:4
    i=ix(j-1); % start with prev pt, find extreme in corresp directio
    while p(i,:)*v(:,j) <= p(i+1,:)*v(:,j)
        i=i+1;
    end
    ix(j) = i;
end

term1 = ix(3);
term = n+ix(1);
ox=ix;
while ix(1)<=term
    ln=norm(u);
    if ln==0
        disp('duplicate point');
    end
    u = u/ln;
    v = [-perp(u);u;perp(u);-u]'; % directions for extremes
    d = [p(ix(1),:)*v(:,1)+p(ix(3),:)*v(:,3)
        p(ix(2),:)*v(:,2)+p(ix(4),:)*v(:,4)];
    [sq,xy] = max(d);
    if result-sq >= -1e-6
        if result-sq > 1e-6
            clf
            hold on
            plot(p(:,1),p(:,2),'ko-'); axis equal
        end
        result = sq;
        fourpts = false;
        delta = abs(diff(d))/2;
        q = [p(ix(1),:); p(ix(2),:); p(ix(3),:); p(ix(4),:)];
        idx=[3-xy;5-xy];
        q(idx,:)=q(idx,:)+delta*v(:,idx)';
        idx=[2:4 1];
        q=q+(repmat(sum((q(idx,:)-q)'.*v(:,idx)),2,1).*v(:,idx))';
        q(end+1,:)=q(1,:);
        plot(q(:,1),q(:,2),'r:d');
    end
    
    ac = p(ix(1),:)-p(ix(3),:);% check 4-point square construction
    u = p(ix(2),:)+ [ac(2) -ac(1)]-p(ix(4),:);
    len = norm(u);
    if len>0
        u = u/len;
        [sq1,i1]=min(-p(:,1)*u(2)+p(:,2)*u(1));
        [sq3,i3]=max(-p(:,1)*u(2)+p(:,2)*u(1));
        [sq2,i2]=max(p(:,1)*u(1)+p(:,2)*u(2));
        [sq4,i4]=min(p(:,1)*u(1)+p(:,2)*u(2));
        [sq,xy] = max([sq3-sq1;sq2-sq4]);
        if result-sq >= -1e-6
            if result-sq > 1e-6
                clf
                hold on
                plot(p(:,1),p(:,2),'ko-'); axis equal
            end
            v = [-perp(u);u;perp(u);-u]'; % directions for extremes
            
            result = sq;
            fourpts = false;
            q = [p(i1,:); p(i2,:); p(i3,:); p(i4,:)];
            idx=[2:4 1];
            q=q+(repmat(sum((q(idx,:)-q)'.*v(:,idx)),2,1).*v(:,idx))';
            q(end+1,:)=q(1,:);
            plot(q(:,1),q(:,2),'m-.d');
        end
    end
    k= 0;
    u = [0,0];
    for j = 1:4 % find next calipers direction
        if (p(ix(j)+1,:)-p(ix(j),:))*(u') >= 0
            k=j;
            u = p(ix(k)+1,:)-p(ix(k),:);
        else
            u=perp(u);
        end
    end
    u = perp(u);
    ix(k)=ix(k)+1;
end

end

