function matrixout2=binMatrix(A,binningsize)
if binningsize==1
    matrixout2=A;
else

    matrixSize=size(A);
    finalrow=floor(matrixSize(1)/binningsize);
    finalcol=floor(matrixSize(2)/binningsize);

    if(mod(matrixSize(1),binningsize)~=0)
        finalrow=finalrow+1;
    end

    if(mod(matrixSize(2),binningsize)~=0)
        finalcol=finalcol+1;
    end

    matrixout1=zeros(finalrow,matrixSize(2));
    matrixout2=zeros(finalcol,finalrow);
    %contract rows
    for i=1:finalrow    
        if i*binningsize<=matrixSize(1)
            matrixout1(i,:)=sum(A((i-1)*binningsize+1:i*binningsize,:));
        else
            if(mod(matrixSize(1),binningsize)>1)
                matrixout1(i,:)=sum(A((i-1)*binningsize+1:end,:));
            else
                matrixout1(i,:)=A(end,:);
            end
        end
    end
    matrixout1=matrixout1';
    %contract cols
    for i=1:finalcol    
        if i*binningsize<=matrixSize(2)
            matrixout2(i,:)=sum(matrixout1((i-1)*binningsize+1:i*binningsize,:));
        else
            if(mod(matrixSize(2),binningsize)>1)
                matrixout2(i,:)=sum(matrixout1((i-1)*binningsize+1:end,:));
            else
                matrixout2(i,:)=matrixout1(end,:);
            end        
            %matrixout2(i,:)=sum(matrixout1((i-1)*binningsize+1:end,:));
        end
    end

    %matrixout2=matrixout2';
    matrixout2=matrixout2';
end