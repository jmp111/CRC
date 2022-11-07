function [U,D,V,X,vari,F]=recursivePCA(X,nc,varargin)
    
    if nargin==0
        disp('####################################################')
        disp('#                  Recursive PCA                   #')
        disp('####################################################')
        disp('# INFO                                             #')
        disp('#   NIPALS singular value decomposition: X=U*D*V''  #')
        disp('#                                                  #')
        disp('# USAGE                                            #')
        disp('#   [U,D,V,Xa,v]=recursivePCA(X,c);                #')
        disp('#   [U,D,V]=recursivePCA(X,c);                     #')
        disp('#   [U,D,V]=recursivePCA(X);                       #')
        disp('#                                                  #')
        disp('# INPUT                                            #')
        disp('#   X:   centered X matrix                   [n,p] #')
        disp('#   c:   number of components to determine   [1,1] #')
        disp('#                                                  #')
        disp('# OUTPUT                                           #')
        disp('#   U:   scores                              [n,c] #')
        disp('#   D:   singular values                     [c,c] #')
        disp('#   V:   loadings                            [p,c] #')
        disp('#   Xa:  X after removal of c components     [n,p] #')
        disp('#   v:   percentage variance explained       [c,1] #')
        disp('#   F:   %variance per feature per component [p,c] #')
        disp('####################################################')
        disp('#                J.M.P. 19/12/2011                 #')
        disp('####################################################')
        return
    end
    
    if nargin<2 || isempty(nc)
        nc=Inf;
    end
    if size(varargin,2)==2
        Elim=varargin{1};
        vari=varargin{2};
        how=2;
    elseif size(varargin,2)==1
        Elim=varargin{1};
        how=1;
    else
        Elim=1e-10;
        how=1;
    end
    [n,p]=size(X);
    nc=min([n,p,nc]);
    SSX=sum(sum(X.^2));
    
    U=zeros(n,nc);
    V=zeros(p,nc);
    F=zeros(p,nc);
    D=zeros(nc,1);
    
    if how==1
        for k=1:nc
            E=Inf;
            v0=zeros(p,1);
            %More elegant, parity check needed?:
    %         vX=var(X);
    %         u=X(:,find(vX==max(vX),1));
            u=X(:,1);
            while E>Elim
                v=(X'*u)/(u'*u);
                v=v/sqrt(v'*v);
                u=(X*v)/(v'*v);
                E=(v0-v);
                v0=v;
                E=E'*E;
            end
            Xp=u*v';
            X=X-Xp;
            F(:,k)=sum(Xp.^2);
            D(k)=sum(F(:,k));
            F(:,k)=F(:,k)/D(k);
            U(:,k)=u/sqrt(D(k));
            V(:,k)=v;
        end
    elseif how==2
        k=1;
        while (sum(D)/SSX)<(vari/100) && k<nc
            E=Inf;
            v0=zeros(p,1);
            %More elegant, parity check needed?:
    %         vX=var(X);
    %         u=X(:,find(vX==max(vX),1));
            u=X(:,1);
            while E>Elim
                v=(X'*u)/(u'*u);
                v=v/sqrt(v'*v);
                u=(X*v)/(v'*v);
                E=(v0-v);
                v0=v;
                E=E'*E;
            end
            Xp=u*v';
            X=X-Xp;
            F(:,k)=sum(Xp.^2);
            D(k)=sum(F(:,k));
            F(:,k)=F(:,k)/D(k);
            U(:,k)=u/sqrt(D(k));
            V(:,k)=v;
            k=k+1;
        end
        V=V(:,1:(k-1));
        D=D(1:(k-1));
        U=U(:,1:(k-1));
        F=F(:,1:(k-1));
    end
    
    vari=D/SSX*100;
    F=F*100;
    D=diag(sqrt(D));
    
end