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
    sq = max(p(ix(1),:)*v(:,1)+p(ix(3),:)*v(:,3),...
        p(ix(2),:)*v(:,2)+p(ix(4),:)*v(:,4));
    if result > sq
        if ix(1)>term1
              disp([ox; ix]);
                disp([result sq]);
            plot(p(:,1),p(:,2),'o-'); axis equal
            disp('wrap');
        end
        result = sq;
        fourpts = false;
    end
    
    ac = p(ix(1),:)-p(ix(3),:);% check 4-point square construction
    u = p(ix(2),:)+ [ac(2) -ac(1)]-p(ix(4),:);
    len = norm(u);
    if len>0
        u = u/len;
        sq = max(max(-p(:,1)*u(2)+p(:,2)*u(1))-min(-p(:,1)*u(2)+p(:,2)*u(1)),...
            max(p(:,1)*u(1)+p(:,2)*u(2))-min(p(:,1)*u(1)+p(:,2)*u(2)));
        if result > sq
            if ix(1)>term1
                disp([ox; ix]);
                disp([result sq]);
                plot(p(:,1),p(:,2),'o-'); axis equal
                disp('wrap');
            end
            result = sq;
            fourpts = true;
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

