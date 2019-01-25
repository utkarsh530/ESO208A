clc;
z =[];
%estimation of eigenvalues
filename = input('Enter the name of the file(with .txt):','s');
fileID = fopen(filename,'r');
line = fgetl(fileID);
n = sscanf(line,'%f');
A = zeros(n ,n);
for i=1:n
    line = fgetl(fileID);
    A(i,1:n) = sscanf(line,'%f');
end
line = fgetl(fileID);
tol = sscanf(line,'%f');
tol = tol/100;
prompt = input('1.Only Largest Eigen Value\n2.All Eigen Value\nEnter your choice(1-2): ');
if(prompt==1)
    err = [];
    i=1;
    z(1:n,1)=zeros(n,1);
    z(1,1)=1;
    y(1:n,1)=z(1:n,1)/norm(z(1:n,1));
    while true
        z(1:n,i+1)=A*y(1:n,i);
        y(1:n,i+1)=z(1:n,i+1)/norm(z(1:n,i+1));
        for j=1:n-1
        err(i)=max(abs((y(j,i+1)-y(j,i))/y(j,i)),abs((y(j+1,i+1)-y(j+1,i))/y(j+1,i)));
        end
        if(err(i)<=tol)
            break;
        end
        i=i+1;
    end
    lambda=(y(1:n,i+1)')*A*y(1:n,i+1)/((y(1:n,i+1)')*y(1:n,i+1));
    fileoID = fopen('output.txt','w');        
    fprintf(fileoID ,'Largest Eigenvalue\n%f\r\n',lambda);
    fclose(fileID);
    fclose(fileoID);
end
if(prompt==2)
    c=1;
    for k=1:n
        eigen(k,c)=A(k,k);
    end    
    while true
        z(1:n,1)=A(1:n,1);
        y(1:n,1)=z(1:n,1)/norm(z(1:n,1));
        for k =1:n-1
            s=zeros(n,1);
            for i=1:k
                s = s + (A(1:n,k+1)'*y(1:n,i))*y(1:n,i);
            end    
            z(1:n,k+1)=A(1:n,k+1)-s;
            y(1:n,k+1)=z(1:n,k+1)/norm(z(1:n,k+1));
        end
        r = zeros(n);
        for i=1:n
            for j=1:n
                if(i<=j)
                    r(i,j)=y(1:n,i)'*A(1:n,j);
                end    
            end    
        end
        A=r*y;
        c=c+1;
        for k=1:n
            eigen(k,c)=A(k,k);
        end
        for k=1:n-1
            err(c-1)=max(abs((eigen(k,c)-eigen(k,c-1))/eigen(k,c-1)),abs((eigen(k+1,c)-eigen(k+1,c-1))/eigen(k+1,c-1)));
        end
        if(err(c-1)<=tol)
            break;
        end    
    end
    fileoID = fopen('output.txt','w');        
    fprintf(fileoID ,'Eigenvalues\n');
    fprintf(fileoID,'%f ',eigen(:,c));
    fclose(fileID);
    fclose(fileoID);    
end    