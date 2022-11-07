function pfdr=bhfdr(p)
%% Benjamini-Hochberg False Discovery Rate
% Y Benjamini, Y Hochberg, (1995) J Stat Roy Soc B Met 57(1):289-300
% pfdr are the FDR adjusted p-values
% JMP - forgot what date, somewhere in August 2012
    
    if ~isvector(p)
        [a,b]=size(p);
        p=p(:);
    else
        a=[];b=[];
    end
    
    n=length(p);
    [p,ii]=sort(p);
    sf=zeros(size(p));
    sf(1:end)=1:n; % faster than sort(ii)
    pfdr=(p*n)./sf;
    
    for k=(n-1):-1:1 % loops back
%         fdr(k)=min(fdr(k:(k+1))); %slower
        if pfdr(k)>pfdr(k+1)
            pfdr(k)=pfdr(k+1);
        end
    end
    
    % while loop is also slower
%     k=1;
%     while k==1
%         up=find(diff(fdr)<0);
%         if ~isempty(up)
%             fdr(up)=fdr(up+1);
%         else
%             k=0;
%         end
%     end
    
    pfdr(ii)=pfdr; % sort like p
    
    if ~isempty(a) && ~isempty(b)
        pfdr=reshape(pfdr,a,b);
    end
    
end