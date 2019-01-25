clc;
%Problem 9 of PS 7
h = input('Grid size h?: ');
% Taking n
n = (2/h)+1;
%Initialising Rearanged func. for Thomas Algo
a = @(x) 1-((-(x+3)/(x+1))*h)/2;
b = @(x) (h^2)*((x+3)/(x+1)^2)-2;
c = @(x) 1+((-(x+3)/(x+1))*h)/2;
f  = @(x) (2*(x+1)+3*((x+3)/(x+1)^2))*(h^2);
prompt = input('1.2nd order Backward Difference\n2.2nd order central difference with Ghost Node\nEnter Your Option(1-2): ');
if(prompt == 2)
        alpha = zeros(1,n-1);
        l = [];
        %Setting l,d,u,B for Thomas
        for i = 1:n-1
            if(i==1)
                l(i) = 0;    
            elseif(i==n-1)
                l(i) = 2;
            else
                l(i) = a(h*(i));
            end    
        end
        u = [];
        for i = 1:n-1
            if(i==n-1)
                u(i) = 0;
            else
                u(i) = c(h*(i));
            end    
        end
        for i = 1:n-1
            d(i) = b(h*(i));    
        end
        B = [];
        for i = 1:n-1
            if(i==1)
                B(i)=f(h*(i))-(a(h*(i)))*5;
            else
                B(i)=f(h*(i));
            end     
        end
        %Finding Alpha and Beta
        alpha(1) = d(1);
        beta = zeros(1,n-1);
        beta(1) = B(1);
        for i=2:n-1
            alpha(i) = d(i)-(l(i)/alpha(i-1))*u(i-1);
            beta(i) = B(i)-(l(i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,n-1);
        x(n-1)=beta(n-1)/alpha(n-1);
        for i = 1:n-2
            x(n-1-i) = (beta(n-1-i)-u(n-1-i)*x(n-i))/alpha(n-1-i);
        end
        X = [];
        for i =1:n
            if(i==1)
                X(i)=5;
            else    
                X(i) = x(i-1);
            end    
        end
        %Plot
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'Points(T0,T1,..TN)\n');
        fprintf(fileoID,'%f ',X);
        fclose(fileoID);
        scatter(0:h:2,X,'filled');
        title('Scatter Plot');
        xlabel('h');
        ylabel('T');
end        
if(prompt==1)
    %Making Coefficient Matrix
    A = zeros(n,n-1);
    B = zeros(n-1,1);
    for i=1:n-1
        if(i==1)
            A(i,i) = b(h*(i));
            A(i,i+1) = c(h*(i));
            A(n,i) = f(h*(i))-(a(h*(i)))*5;
        elseif(i==n-1)
            A(i,i-2) = 1/3;
            A(i,i-1) = -4/3;
            A(i,i) = 1;
            A(n,i)  = 0;          
        else
            A(i,i-1) = a(h*i);
            A(i,i) = b(h*i);
            A(i,i+1) = c(h*i);
            A(n,i) = f(h*i);            
        end
    end
    n = n-1;
    % Applying Gauss Jordan with Pivoting
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
        X = [];
        for i =1:n+1
            if(i==1)
                X(i)=5;
            else    
                X(i) = x(i-1);
            end    
        end
        %Plot
        fileoID = fopen('output.txt','w');
        fprintf(fileoID,'Points(T0,T1,..TN)\n');
        fprintf(fileoID,'%f ',X);
        fclose(fileoID);
        scatter(0:h:2,X,'filled');
        title('Scatter Plot');
        xlabel('h');
        ylabel('T');
end    