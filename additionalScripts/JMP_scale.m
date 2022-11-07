function [tr,te,varargout]=JMP_scale(tr,te,scaling,varargin)
% [tr,te,varargout]=...
% [tr,te,mc,va]=... % for 'r' or 'c'
% [tr,te,mc,mc1,mc2]=... % for 'b', where mc1 for columns, mc2 for rows
    
    [ntr,ptr]=size(tr);
    [nte,pte]=size(te);
    scond=0;
    if ptr~=pte && ~isempty(te)
        scond=1;
    end
    
    if scond==0 && ischar(scaling)
        if strcmp(scaling,'mc')
            lambda=0;
        elseif strcmp(scaling,'pa')
            lambda=0.5;
        elseif strcmp(scaling,'as')
            lambda=1;
        else
            scond=1;
        end
    elseif scond==0 && isnumeric(scaling)
        if numel(scaling)==1
            lambda=max([0,min([1,scaling])]);
        else
            scond=1;
        end
    else
        scond=1;
    end
    
    if scond==0 && nargin==4
        direction=varargin{1};
        if ischar(direction) && (strcmp(direction(1),'c') || strcmp(direction(1),'r') || strcmp(direction(1),'b'))
            direction=direction(1);
        else
            scond=1;
        end
    elseif scond==0
        direction='c';
    end
    
    if scond==0 && strcmp(direction,'b') && numel(lambda)~=1
        scond=1;
    end
    
    if scond==0
        if strcmp(direction,'c')
            mc=mean(tr,1);
            tr=tr-mc(ones(ntr,1),:);
            if ~isempty(te)
                te=te-mc(ones(nte,1),:);
            end
            va=std(tr,0,1).^lambda;
            va(va==0)=1;
            tr=tr./va(ones(ntr,1),:);
            if ~isempty(te)
                te=te./va(ones(nte,1),:);
            end
        elseif strcmp(direction,'r')
            mc=mean(tr,2);
            tr=tr-mc(:,ones(1,ptr));
            if ~isempty(te)
                te=te-mc(:,ones(1,pte));
            end
            va=std(tr,0,2).^lambda;
            va(va==0)=1;
            tr=tr./va(:,ones(1,ptr));
            if ~isempty(te)
                te=te./va(:,ones(1,pte));
            end
        elseif strcmp(direction,'b')
            mc1=mean(tr,1);
            mc2=mean(tr,2);
            mc=mean(tr(:));
            tr=tr-mc1(ones(ntr,1),:)-mc2(:,ones(1,ptr))+mc;
            if ~isempty(te)
                te=te-mc1(ones(nte,1),:)-mc2(:,ones(1,pte))+mc;
            end
        end
    end
    
    for k=3:nargout
        if scond==1
            if isodd(k)
                varargout{k-2}=0;
            else
                varargout{k-2}=1;
            end
        else
            if strcmp(direction,'r') || strcmp(direction,'c')
                if k<=4
                    if isodd(k)
                        varargout{k-2}=mc;
                    else
                        varargout{k-2}=va;
                    end
                else
                    if isodd(k)
                        varargout{k-2}=0;
                    else
                        varargout{k-2}=1;
                    end
                end
            else
                if k<=5
                    if k==3
                        varargout{k-2}=mc;
                    elseif k==4
                        varargout{k-2}=mc1;
                    elseif k==5
                        varargout{k-2}=mc2;
                    end
                else
                    if isodd(k)
                        varargout{k-2}=0;
                    else
                        varargout{k-2}=1;
                    end
                end
            end
        end
    end
    
end