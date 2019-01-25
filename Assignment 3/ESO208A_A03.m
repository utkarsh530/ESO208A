clc;
filename = input('Enter the name of the file(with .txt):','s');
fileID = fopen(filename,'r');
formatSpec = '%f %f';
sizeI = [2 Inf];
I = fscanf(fileID,formatSpec,sizeI);
fclose(fileID);
[n,m]=size(I);
disp('What you want to do?')
prompt = input('Matrix Operations\n1. Fit a least square polynomial\n2. Fit a Lagrange Interpolation Polynomial\n3. Fit Cubic Splines\nEnter your option(1-3):');
if(prompt == 1)
    order = input('Order of the polynomial (<n)?: \n');
    intercept = input('Should Intercept be zero?(Y/N) \n', 's');
    A = zeros(order+1,order+1);
    for i=1:order+1
        for j=1:order+1     
            for k=1:m
                A(i,j) = A(i,j) + ((I(1,k))^(i-1) * (I(1,k))^(j-1));
            end
        end
    end
    B = zeros(order+1,1);
    for i=1:order+1
        for j=1:m
            B(i,1) = B(i,1) + I(2,j)*(I(1,j)^(i-1));
        end
    end
    C = A\B;
    C = flipud(C);
    P = poly2sym(C);
    scatter(I(1,:) , I(2,:))
    grid on
    title('Scatter Plot');
    hold on
    ezplot(P,[min(I(1,:))-10,max(I(1,:))+10])
    title('Curve');
end
if(prompt == 2)
    globalcoeff = zeros(m,1);
    for i=1:m
        it = 1;
        deno = 1;
        root = zeros(m-1,1);
        for j=1:m
            if j~=i
                root(it,1) = I(1,j);
                it = it + 1;
                deno = deno*(I(1,i)-I(1,j));
            end    
        end
        deno = 1/deno;
        deno = I(2,i)*deno;
        p = poly(root);
        for k=1:m
            globalcoeff(k,1) = globalcoeff(k,1) + p(1,k)*deno; 
        end    
    end
    disp('Polynomial is: ');
    disp(poly2sym(globalcoeff));
    scatter(I(1,:) , I(2    ,:));
    grid on
    axis([-inf  inf    -inf  inf])
    hold on
    ezplot(poly2sym(globalcoeff));
    axis([-inf  inf    -inf  inf])
    title('Curve and Scatter Plot');    
end  
if(prompt == 3)
    opt = input('1) for Natural Spline\n 2) for Not-a-knot\n 3) for Periodic\n 4) Clamped Spline \nEnter(1-4): ');
    if(opt==1)
        v = zeros(m,1);
        h = [];
        for i = 1:m-1
            h(i) = I(1,i+1)-I(1,i);
        end
        alpha = zeros(1,m-2);
        l = [];
        for i = 1:m-2
            if(i==1)
                l(i) = 0;
            else
                l(i) = h(i);
            end    
        end
        u = [];
        for i = 2:m-1
            if(i==m-1)
                u(i-1) = 0;
            else
                u(i-1) = h(i);
            end    
        end
        for i = 1:m-2
            d(i) = 2*(h(i)+h(i+1));    
        end
        b = [];
        for i = 1:m-2
            b(i)=3*((I(2,i+2)-I(2,i+1))/(I(1,i+2)-I(1,i+1))-((I(2,i+1)-I(2,i)))/(I(1,i+2)-I(1,i+1)));
        end    
        alpha(1) = d(1);
        beta = zeros(1,m-2);
        beta(1) = b(1);
        for i=2:m-2
            alpha(i) = d(i)-(l(i)/alpha(i-1))*u(i-1);
            beta(i) = b(i)-(l(i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,m-2);
        x(m-2)=beta(m-2)/alpha(m-2);
        for i = 1:m-3
            x(m-2-i) = (beta(m-2-i)-u(m-2-i)*x(m-1-i))/alpha(m-2-i);
        end
        v = [];
        for i=1:m
            if(i==1 || i==m)
            v(i)=0;
            else
                v(i)=x(i-1);
            end
        end    
        U = [];
        for i = 1:m-1
            U(i)=((I(2,i+1)-I(2,i))/(I(1,i+1)-I(1,i)))-(h(i)*(v(i+1)+2*v(i)))/3;
        end
        U(m) = 3*((I(2,m)-I(2,m-1))/(I(1,m)-I(1,m-1)))-2*U(m-1)-(h(m-1)*v(m-1));
        p =[];
        coeffmat = [];
        syms f(x);
        fileoID = fopen('output.txt','w');        
        for i=1:m-1
            f(x) =I(2,i)+(x-I(1,i))*U(i)+((x-I(1,i))^2)*v(i)+((x-I(1,i))^3)*(v(i+1)-v(i))/(3*h(i));        
            fprintf(fileoID ,'\nPolynomial\n\r\n');
            fprintf(fileoID,'%s  ',f(x));
            fprintf(fileoID,'\r\n');            
            hold on
            ezplot(f(x),[I(1,i),I(1,i+1)])
            hold on            
        end
        grid on
        axis([-inf  inf    -inf  inf]);
        scatter(I(1,:) , I(2    ,:));
    end  
    if(opt==2)

        v = zeros(m,1);
        h = [];
        for i = 1:m-1
            h(i) = I(1,i+1)-I(1,i);
        end
        alpha = zeros(1,m-2);
        l = [];
        for i = 1:m-2
            if(i==1)
                l(i) = 0;
            else
                l(i) = h(i);
            end    
        end
        u = [];
        for i = 2:m-1
            if(i==m-1)
                u(i-1) = 0;
            else
                u(i-1) = h(i);
            end    
        end
        for i = 1:m-2
            d(i) = 2*(h(i)+h(i+1));    
        end
        b = [];
        for i = 1:m-2
            b(i)=3*((I(2,i+2)-I(2,i+1))/(I(1,i+2)-I(1,i+1))-((I(2,i+1)-I(2,i)))/(I(1,i+2)-I(1,i+1)));
        end    
        alpha(1) = d(1);
        beta = zeros(1,m-2);
        beta(1) = b(1);
        for i=2:m-2
            alpha(i) = d(i)-(l(i)/alpha(i-1))*u(i-1);
            beta(i) = b(i)-(l(i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,m-2);
        x(m-2)=beta(m-2)/alpha(m-2);
        for i = 1:m-3
            x(m-2-i) = (beta(m-2-i)-u(m-2-i)*x(m-1-i))/alpha(m-2-i);
        end
        v = [];
        for i=1:m
            if(i==1)
            v(i)=x(1)-(h(1)/h(2))*(x(2)-x(1));
            end
            if(i==m)
            v(i)= x(m-2)+(h(m-1)/h(m-2))*(x(m-2)-x(m-3));
            end
            if (i~=1 && i~=m)
                v(i)=x(i-1);
            end
        end    
        U = [];
        for i = 1:m-1
            U(i)=((I(2,i+1)-I(2,i))/(I(1,i+1)-I(1,i)))-(h(i)*(v(i+1)+2*v(i)))/3;
        end
        U(m) =3*((I(2,m)-I(2,m-1))/(I(1,m)-I(1,m-1)))-2*U(m-1)-(h(m-1)*v(m-1));
        p =[];
        syms f(x);
        fileoID = fopen('output.txt','w');
        for i=1:m-1
            f(x) =I(2,i)+(x-I(1,i))*U(i)+((x-I(1,i))^2)*v(i)+((x-I(1,i))^3)*(v(i+1)-v(i))/(3*h(i));        
            fprintf(fileoID ,'\nPolynomial\n\r\n');
            fprintf(fileoID,'%s  ',f(x));
            fprintf(fileoID,'\r\n');            
            hold on
            ezplot(f(x),[I(1,i),I(1,i+1)])
            hold on
        end
        grid on
        axis([-inf  inf    -inf  inf]);
        scatter(I(1,:) , I(2    ,:));
        fclose(fileoID);
    end
    if(opt==3)
        v = zeros(m,1);
        h = [];
        for i = 1:m-1
            h(i) = I(1,i+1)-I(1,i);
        end
        alpha = zeros(1,m-2);
        l = [];
        for i = 1:m-2
            if(i==1)
                l(i) = 0;
            else
                l(i) = h(i);
            end    
        end
        u = [];
        for i = 2:m-1
            if(i==m-1)
                u(i-1) = 0;
            else
                u(i-1) = h(i);
            end    
        end
        for i = 1:m-2
            d(i) = 2*(h(i)+h(i+1));    
        end
        b = [];
        for i = 1:m-2
            b(i)=3*((I(2,i+2)-I(2,i+1))/(I(1,i+2)-I(1,i+1))-((I(2,i+1)-I(2,i)))/(I(1,i+2)-I(1,i+1)));
        end    
        alpha(1) = d(1);
        beta = zeros(1,m-2);
        beta(1) = b(1);
        for i=2:m-2
            alpha(i) = d(i)-(l(i)/alpha(i-1))*u(i-1);
            beta(i) = b(i)-(l(i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,m-2);
        x(m-2)=beta(m-2)/alpha(m-2);
        for i = 1:m-3
            x(m-2-i) = (beta(m-2-i)-u(m-2-i)*x(m-1-i))/alpha(m-2-i);
        end
        v = [];
        for i=1:m
            if(i==1)
            v(i)=1;
            end
            if(i==m)
            v(i)= 1;
            end
            if (i~=1 && i~=m)
                v(i)=x(i-1);
            end
        end    
        U = [];
        for i = 1:m-1
            U(i)=((I(2,i+1)-I(2,i))/(I(1,i+1)-I(1,i)))-(h(i)*(v(i+1)+2*v(i)))/3;
        end
        U(m-1) = 3*((I(2,m-1)-I(2,m-2))/(I(1,m-1)-I(1,m-2)))-(2*U(m-2))-(h(m-2)*v(m-2));
        U(m) = 3*((I(2,m)-I(2,m-1))/(I(1,m)-I(1,m-1)))-(2*U(m-1))-(h(m-1)*v(m-1));
        U(1) = U(m);
        v(1) = (((3/h(1))*(((I(2,2)-I(2,1))/(I(1,2)-I(1,1)))-U(1)))-v(2))/2;
        v(m) = v(1);
        p =[];
        syms f(x);
        fileoID = fopen('output.txt','w');
        for i=1:m-1
            f(x) =I(2,i)+(x-I(1,i))*U(i)+((x-I(1,i))^2)*v(i)+((x-I(1,i))^3)*(v(i+1)-v(i))/(3*h(i));        
            fprintf(fileoID ,'\nPolynomial\n\r\n');
            fprintf(fileoID,'%s  ',f(x));
            fprintf(fileoID,'\r\n');            
            hold on
            ezplot(f(x),[I(1,i),I(1,i+1)])
            hold on
        end
        grid on
        axis([-inf  inf    -inf  inf]);
        scatter(I(1,:) , I(2    ,:));
        fclose(fileoID);        
    end
    if(opt==4)
        fprintf('enter first derivative value \n');
        s1=input('');
        fprintf('enter second derivative value \n');
        s2=input('');        
        v = zeros(m,1);
        h = [];
        for i = 1:m-1
            h(i) = I(1,i+1)-I(1,i);
        end
        alpha = zeros(1,m-2);
        l = [];
        for i = 1:m-2
            if(i==1)
                l(i) = 0;
            else
                l(i) = h(i);
            end    
        end
        u = [];
        for i = 2:m-1
            if(i==m-1)
                u(i-1) = 0;
            else
                u(i-1) = h(i);
            end    
        end
        for i = 1:m-2
            d(i) = 2*(h(i)+h(i+1));    
        end
        b = [];
        for i = 1:m-2
            b(i)=3*((I(2,i+2)-I(2,i+1))/(I(1,i+2)-I(1,i+1))-((I(2,i+1)-I(2,i)))/(I(1,i+2)-I(1,i+1)));
        end    
        alpha(1) = d(1);
        beta = zeros(1,m-2);
        beta(1) = b(1);
        for i=2:m-2
            alpha(i) = d(i)-(l(i)/alpha(i-1))*u(i-1);
            beta(i) = b(i)-(l(i)/alpha(i-1))*beta(i-1);
        end
        x = zeros(1,m-2);
        x(m-2)=beta(m-2)/alpha(m-2);
        for i = 1:m-3
            x(m-2-i) = (beta(m-2-i)-u(m-2-i)*x(m-1-i))/alpha(m-2-i);
        end
        v = [];
        for i=1:m
            if(i==1)
            v(i)=1;
            end
            if(i==m)
            v(i)= 1;
            end
            if (i~=1 && i~=m)
                v(i)=x(i-1);
            end
        end    
        U = [];
        for i = 1:m-1
            U(i)=((I(2,i+1)-I(2,i))/(I(1,i+1)-I(1,i)))-(h(i)*(v(i+1)+2*v(i)))/3;
        end
        U(1)=s1;
        U(m)=s2;
        U(m-1)=(U(m) - 3*((I(2,m)-I(2,m-1))/(I(1,m)-I(1,m-1)))+(h(m-1)*v(m-1)))/(-2);
        v(1) = (((3/h(1))*(((I(2,2)-I(2,1))/(I(1,2)-I(1,1)))-U(1)))-v(2))/2;
        v(m)= ((3/h(m-1))*(((I(2,m)-I(2,m-1))/(I(1,m)-I(1,m-1)))-U(m-1)))-2*v(m-1);
        syms f(x);
        fileoID = fopen('output.txt','w');
        for i=1:m-1
            f(x) =I(2,i)+(x-I(1,i))*U(i)+((x-I(1,i))^2)*v(i)+((x-I(1,i))^3)*(v(i+1)-v(i))/(3*h(i));        
            fprintf(fileoID ,'\nPolynomial\n\r\n');
            fprintf(fileoID,'%s  ',f(x));
            fprintf(fileoID,'\r\n');            
            hold on
            ezplot(f(x),[I(1,i),I(1,i+1)])
            hold on
        end
        grid on
        axis([-inf  inf    -inf  inf]);
        scatter(I(1,:) , I(2    ,:));
        fclose(fileoID);         
    end    
end    