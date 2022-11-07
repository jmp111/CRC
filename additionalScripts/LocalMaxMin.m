function lm=LocalMaxMin(v)
    
    n=size(v,1); %if n=1 then lm needs to be transposed
    v=v(:);
    lm=zeros(size(v));
    v2=[0;v;0];
    d=diff(v2);
    d=[d(1:(end-1)) d(2:end)];
    s=sign(d);
    s=s(:,1).*s(:,2);
    s=(s==-1);
    lm(s)=1;
    lm=-double(d(:,1)<d(:,2)).*lm+double(d(:,1)>d(:,2)).*lm;
    if n==1;
        lm=lm';
    end
    
end