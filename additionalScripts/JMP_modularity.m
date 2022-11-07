function [m,th,mmb,mb]=JMP_modularity(D,perc,m2,p)
    
    D2=D(p,p);
    D2(D2<perc)=0;
    t1=tril(ones(size(D)),-1);
    D2=D2.*(t1+t1'); % sets diagonal to zero
%     D2(D2~=0)=1; % conventional modularity 0/1, presence/absence of edge. Comment this line if the size of edge should be taken into account
    
    m=Modul(D2,m2,p);
    
    if nargout>=3 % question is, can I use m2, or should I resample the indices?
        n=1000;
        b=D2(t1==1);
        v=length(b);
        mb=zeros(n,size(D2,1));
        parfor k=1:n %Bootstrap estimate of modularity
            s=randsample(b,v,'true');
            S=squareform(s);
            mb(k,:)=Modul(S,m2,sort(p));
        end
        mmb=max(mb);
        thboot=m.*double(m>mmb);
        th=imax(thboot);
    elseif nargout==2
        th=imax(m);
    end
     
end

function m=Modul(D2,m2,p)
    
    m=zeros(1,size(D2,1));
    lbo=sum(sum(tril(D2))); %number of significant links between objects
    for k=1:size(D2,1)
        bla=cluster(m2,'maxclust',k);
        bla=bla(p);
        for kk=1:k
            cmat=D2(bla==kk,bla==kk);
            eins=sum(sum(tril(cmat)));%number of edges in cluster kk
            kkdegree=sum(sum(cmat));%sum of degrees of all objects in cluster kk
            m(k)=m(k)+((eins/lbo)-(kkdegree/(2*lbo))^2);
        end
    end
    
end