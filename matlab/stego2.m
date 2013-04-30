function res = stego2(X)
% for fminsearch to find the stego that needs the largest bounding square
% to cover.
%[X,sq,flag,out]= fminsearch(@stego2, [.15,.6]',optimset('TolX',1e-8))
%[X,sq,flag,out]= fminsearch(@stego2, X,optimset('TolX',1e-8))
% th=X(1); c = X(2)*cos(th)+(1-X(2))*sin(th);
% p = [-c*tan(th) c;
% -sin(th) 2*c-cos(th);
% 0,0; % origin
% cos(th)-sin(th) 2*c-cos(th)-sin(th);
% cos(th) sin(th);
% cos(th)+(sin(th)-c)*tan(th) c];
% x=p(:,1)'; y=p(:,2)';
% plot(x([6 1 3 5 6 4 2 1]'),y([6 1 3 5 6 4 2 1]'),'r:d')
% axis equal

% x(1) = th, x(2) = c
th=X(1); c = X(2)*cos(th)+(1-X(2))*sin(th);
% function res = stego(th,Y)
% c = Y*cos(th)+(1-Y)*sin(th);
p = [-c*tan(th) c;
    -sin(th) 2*c-cos(th); 
    0,0; % origin
    cos(th)-sin(th) 2*c-cos(th)-sin(th);
    cos(th) sin(th);
    cos(th)+(sin(th)-c)*tan(th) c];

res = -enclsquaresize(p);
%michael's coordinates
%disp([1-p(2,:)*[cos(th);sin(th)],1-p(2,:)*[-sin(th);cos(th)],  res])
    
    
    
    