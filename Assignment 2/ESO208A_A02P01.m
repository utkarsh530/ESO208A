clc;
%Solution of linear System
prompt = input('Matrix Operations\n1. Solve a system of Equation\n2. Perform a LU decomposition\n3. Perform a Matrix Inversion\nEnter your option(1-3):');
if(prompt == 1)
    tridia = input('Is the system Tri-diagonal?(Y/N):','s');
    if(tridia == 'Y')
        filename = input('Enter the name of the file(with .txt):','s');
        fileID = fopen(filename,'r');
        line = fgetl(fileID);
        n = sscanf(line,'%f');
        A = zeros(4,n);
        for i=1:4
            line = fgetl(fileID);
            A(i,1:n) = sscanf(line,'%f');
        end
        alpha = zeros(1,n);
        alpha(1) = A(2,1);
        beta = zeros(1,n);
        beta(1) = A(4,1);
        for i=2:n
            alpha(i) = A(2,i)-(A(1,i)/alpha(i-1))*A(3,i-1);
            beta(i) = A(4,i)-(A(1,i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,n);
        x(n)=beta(n)/alpha(n);
        for i = 1:n-1
            x(n-i) = (beta(n-i)-A(3,n-i)*x(n+1-i))/alpha(n-i);
        end
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'Roots(x1,x2,x3,x4..xn)\n');
        fprintf(fileoID,'%f ',x);
        fclose(fileID);
        fclose(fileoID);
    else
        filename = input('Enter the name of the file(with .txt):','s');
        fileID = fopen(filename,'r');
        line = fgetl(fileID);
        n = sscanf(line,'%f');
        A = zeros(n+1,n);
        for i=1:n+1
            line = fgetl(fileID);
            A(i,1:n) = sscanf(line,'%f');
        end
        for k=1:n-1
            for i=k+1:n
                big = A(k,k);
                rswap = k;
                if(big<=A(i,k))
                    big = A(i,k);
                    rswap = i;
                end    
            end
            temp = A(i,:);
            A(i,:) = A(rswap,:);
            A(rswap,:) = temp;
            temp = A(n+1,i);
            A(n+1,i)=A(n+1,rswap);
            A(n+1,rswap) = temp;  
            for i = k+1:n
                A(n+1,i)= A(n+1,i)-(A(i,k)/A(k,k))*A(n+1,k);                
                A(i,:) = A(i,:)- (A(i,k)/A(k,k))*A(k,:);
            end    
        end
        x = zeros(1,n);
        x(n) = A(n+1,n)/A(n,n);
        for i=1:n-1
            s=0;
            for j = n+1-i:n
                s = s + A(n-i,j)*x(j);
            end
            x(n-i) = (A(n+1,n-i)-s)/A(n-i,n-i);
        end
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'Roots(x1,x2,..xn)\n');
        fprintf(fileoID,'%f ',x);
        fclose(fileID);
        fclose(fileoID);
    end  
end
if(prompt == 2)
    cLU = input('Is the matrix is symmetric and positive definite?(Y/N)','s');
    if(cLU == 'Y')
        filename = input('Enter the name of the file(with .txt):','s');
        fileID = fopen(filename,'r');
        line = fgetl(fileID);
        n = sscanf(line,'%f');
        A = zeros(n+1,n);
        for i=1:n+1
            line = fgetl(fileID);
            A(i,1:n) = sscanf(line,'%f');
        end
        B=A(n+1,:);
        A=A(1:n,1:n);
        for k = 1:n-1
            [col_max,r_index]=max(abs(A(k:n,k:n)));
            [max_val,c_index]=max(col_max);
            row(k)=r_index(c_index)+k-1;
            col(k)=c_index+k-1;
            A([k,row(k)],:)=A([row(k),k],:);
            A(:,[k,col(k)])=A(:,[col(k),k]);
            indexes=[row,col];
            if(A(k,k)==0)
                break;
            end    
        end
        L = zeros(n);
        for j=1:n
            s=0;
            for k=1:j-1
                s=s+L(j,k)^2;
            end
            L(j,j)=sqrt(A(j,j)-s);
            for i=j+1:n
                s1=0;
                for k=1:j-1
                    s1=s1+L(i,k)*L(j,k);
                end
                L(i,j)=(A(i,j)-s1)/L(j,j);
            end            
        end  
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'L Matrix(Cholesky Decomposition)\n');
        fprintf(fileoID,'   row  \n');
        fprintf(fileoID,'%6.f ',row);
        fprintf(fileoID,'\n   column  \n');
        fprintf(fileoID,'%6.f ',col);        
        fprintf(fileoID ,'\nMatrix\n\r\n');
        for i=1:1:n
        fprintf(fileoID,'%6.4f  ',L(i,:));
        fprintf(fileoID,'\r\n');
        end
        fclose(fileID);
        fclose(fileoID);
    elseif(cLU=='N')
        filename = input('Enter the name of the file(with .txt):','s');
        fileID = fopen(filename,'r');
        line = fgetl(fileID);
        n = sscanf(line,'%f');
        A = zeros(n+1,n);
        for i=1:n+1
            line = fgetl(fileID);
            A(i,1:n) = sscanf(line,'%f');
        end
        B=A(n+1,:);
        A=A(1:n,1:n);
        for k = 1:n-1
            [col_max,r_index]=max(abs(A(k:n,k:n)));
            [max_val,c_index]=max(col_max);
            row(k)=r_index(c_index)+k-1;
            col(k)=c_index+k-1;
            A([k,row(k)],:)=A([row(k),k],:);
            A(:,[k,col(k)])=A(:,[col(k),k]);
            if(A(k,k)==0)
                break;
            end    
        end
        opt=input('1.Doolittle or 2.Crout(1-2): ');
        if(opt==1)
            L=zeros(n);
            U=zeros(n);
            L(1,1)=1;
            L(2,2)=1;
            L(3,3)=1;
            for k=1:n
                for j=k:n
                    s1=0;
                    for m=1:k-1
                        s1=s1+L(k,m)*U(m,j);
                    end
                    U(k,j)=A(k,j)-s1;
                end
                for i=k+1:n
                    s2=0;
                    for m=1:k-1
                        s2=s2+L(i,m)*U(m,k);
                    end
                    L(i,k)=(A(i,k)-s2)/U(k,k);
                end    
            end
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'LU Matrix(Doolittle Decomposition)\n');
        fprintf(fileoID,'   row  \n');
        fprintf(fileoID,'%6.f ',row);
        fprintf(fileoID,'\n   column  \n');
        fprintf(fileoID,'%6.f ',col);                  
        fprintf(fileoID ,'\nL Matrix\n\r\n');
        for i=1:1:n
            fprintf(fileoID,'%6.4f  ',L(i,:));
            fprintf(fileoID,'\r\n');
        end
        fprintf(fileoID ,'\nU Matrix\n\r\n');
        for i=1:1:n
            fprintf(fileoID,'%6.4f  ',U(i,:));
            fprintf(fileoID,'\r\n');
        end
        else
            L=zeros(n);
            U=zeros(n);
            U(1,1)=1;
            U(2,2)=1;
            U(3,3)=1;
            for k=1:n
                for i=k:n
                    s1=0;
                    for m=1:k-1
                        s1=s1+L(i,m)*U(m,k);
                    end
                    L(i,k)=A(i,k)-s1;
                end
                for j=k+1:n
                    s2=0;
                    for m=1:k-1
                        s2=s2+L(k,m)*U(m,j);
                    end
                    U(k,j)=(A(k,j)-s2)/L(k,k);
                end    
            end
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'LU Matrix(Crout Decomposition)\n');
        fprintf(fileoID,'   row  \n');
        fprintf(fileoID,'%6.f ',row);
        fprintf(fileoID,'\n   column  \n');
        fprintf(fileoID,'%6.f ',col);                  
        fprintf(fileoID ,'\nL Matrix\n\r\n');
        for i=1:1:n
            fprintf(fileoID,'%6.4f  ',L(i,:));
            fprintf(fileoID,'\r\n');
        end
        fprintf(fileoID ,'\nU Matrix\n\r\n');
        for i=1:1:n
            fprintf(fileoID,'%6.4f  ',U(i,:));
            fprintf(fileoID,'\r\n');
        end
        fclose(fileID);
        fclose(fileoID);
        end    
    end    
end
if(prompt==3)
    filename = input('Enter the name of the file(with .txt):','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n);
    for i=1:n
        line = fgetl(fileID);
        A(i,1:n) = sscanf(line,'%f');
    end
    I=zeros(n);
    for i =1:n
        I(i,i)=1;
    end 
    for k=1:n
        I(k,:)=I(k,:)/A(k,k);
        A(k,:)=A(k,:)/A(k,k);
        for i =1:n
            if i~=k
                I(i,:)=I(i,:)-(A(i,k)*I(k,:));
                A(i,:)=A(i,:)-(A(i,k)*A(k,:));
            end    
        end
        fileoID = fopen('output.txt','w');        
        fprintf(fileoID ,'\nInverse Matrix\n\r\n');
        for i=1:1:n
            fprintf(fileoID,'%6.4f  ',I(i,:));
            fprintf(fileoID,'\r\n');
        end
    end
    fclose(fileID);
    fclose(fileoID);
end    