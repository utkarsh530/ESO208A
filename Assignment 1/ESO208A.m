clc;
aer = [];
aes = [];
syms x;
%Finding a root of the non-linear equation ESO208A
prompt = input('Find root of the non-linear equation\nIs the Equation a polynomial?(Y/N):','s');
if(prompt == 'N')
   method = input('1. Bisection\n2. False Position\n3. Fixed-Point\n4. Newton-Raphson\n5. Secant\nEnter your option(1-5): ');
   if(method == 1)
       y = input('Enter f(x): ','s');
       y = evalin(symengine, y);
       a = input('Enter a point(1): ');
       b = input('Enter a point(2): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');
       m = [];
       ae = [];
       for i = 1:itno
           m(i)=(a+b)/2;
           if(abs(subs(y,m(i)))<=valf)
               flag = 1;
               break
           end
           if(subs(y,m(i))*subs(y,a)<0)
               b=m(i);
           else
               a=m(i);
           end
           if(i>1)
               ae(i)=abs((m(i)-m(i-1))/m(i-1));
               if(ae(i)<=(approxerror))
                   flag = 2;
                   break;
               end    
           end    
       end
       fprintf('The value of the root is: %s\n',m(i));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       fplot(y);
       grid on;
       figure;
       plot(1:i,ae,'-o');
       grid on;
   end
   if(method == 2)
       y = input('Enter f(x): ','s');
       y = evalin(symengine, y);
       a = input('Enter a point(1): ');
       b = input('Enter a point(2): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');       
       m = [];
       ae = [];
       for i = 1:itno
           m(i)=a-(subs(y,a)*((b-a)/(subs(y,b)-subs(y,a))));
           if(abs(subs(y,m(i)))<=valf)
               flag = 1;
               break
           end
           if(subs(y,m(i))*subs(y,a)<0)
               b=m(i);
           else
               a=m(i);
           end
           if(i>1)
               ae(i)=abs((m(i)-m(i-1))/m(i-1));
               if(ae(i)<=(approxerror))
                   flag = 2;
                   break
               end    
           end    
       end
       fprintf('The value of the root is: %s\n',m(i));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       fplot(y);
       grid on;
       figure;
       plot(1:i,ae,'-o');
       grid on;
   end
   if(method == 3)
       g = input('Enter a g(x) such that f(x) is expressed as x=g(x): ','s');
       g = evalin(symengine, g);
       a = input('Enter a point(1): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');
       m = [];
       ae = [];
       for i = 1:itno
           if(i>1)
            m(i)=subs(g,m(i-1));
           else
               m(i) = a;
           end
           if(abs(subs(y,m(i)))<=valf)
               flag = 1;
               break
           end
           if(i>1)
               ae(i)=abs((m(i)-m(i-1))/m(i-1));
               if(ae(i)<=(approxerror))
                   flag = 2;
                   break
               end    
           end    
       end
       fprintf('The value of the root is: %s\n',m(i));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       fplot(y);
       grid on;
       figure;
       plot(1:i,ae,'-o');
       grid on;
   end
   if(method == 4)
       y = input('Enter f(x): ','s');
       y = evalin(symengine, y);
       dy = input('Enter f''(x): ','s');
       dy = evalin(symengine, dy);
       a = input('Enter a point(1): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');
       m = [];
       ae = [];
       for i = 1:itno
           if(i>1)
            m(i)=m(i-1)-(subs(y,m(i-1))/subs(dy,m(i-1)));
           else
               m(i) = a;
           end
           if(abs(subs(y,m(i)))<=valf)
               flag = 1;
               break
           end
           if(i>1)
               ae(i)=abs((m(i)-m(i-1))/m(i-1));
               if(ae(i)<=(approxerror))
                   flag = 2;
                   break
               end    
           end    
       end
       fprintf('The value of the root is: %s\n',m(i));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       fplot(y);
       grid on;
       figure;
       plot(1:i,ae,'-o');
       grid on;
   end
   if(method == 5)
       y = input('Enter f(x): ','s');
       y = evalin(symengine, y);
       a = input('Enter a point(1): ');
       b = input('Enter a point(2): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');
       m = [];
       ae = [];
       for i = 1:itno
           if(i>2)
            m(i)=m(i-1)-(subs(y,m(i-1))*(m(i-1)-m(i-2))/(subs(y,m(i-1))-subs(y,m(i-2))));
           else
               m(1) = a;
               m(2) = b;
           end
           if(abs(subs(y,m(i)))<=valf)
               flag = 1;
               break
           end
           if(i>1)
               ae(i)=abs((m(i)-m(i-1))/m(i-1));
               if(ae(i)<=(approxerror))
                   flag = 2;
                   break
               end    
           end    
       end
       fprintf('The value of the root is: %s\n',m(i));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       fplot(y);
       grid on;
       figure;
       plot(1:i,ae,'-o');
       grid on;
   end
elseif(prompt == 'Y')
    poly  = [];
    order = input('order of the polynomial: ');
    for i = 1: order+1
        poly(i) = input('input coeff: ');
    end    
    method = input('1. Muller\n2. Bairstow\nEnter your choice(1-2): ');
    if(method == 1)
       a = input('Enter a point(1): ');
       b = input('Enter a point(2): ');
       c = input('Enter a point(3): ');
       approxerror = input('Enter stopping criteria\nRelative approximate error: ');
       valf = input('How close f(x) is to zero: ');
       itno = input('Maximum Iteration number: ');
       m(1) = a;
       m(2) = b;
       m(3) = c;
       ae = [];
       for i=1:itno
           h0 = m(i+1) - m(i);
           h1 = m(i+2) - m(i+1);
           divd1 = (polyval(poly, m(i+1)) - polyval(poly, m(i)))/h0;
           divd2 = (polyval(poly, m(i+2)) - polyval(poly, m(i+1)))/h1;
           a = (divd2 - divd1)/(h0+h1);
           b = a*h1 + divd2;
           c = polyval(poly,m(i+2));
           deno1 = abs(b + sqrt(b^(2) - 4*a*c));
           deno2 = abs(b - sqrt(b^(2) - 4*a*c));
           if(deno1 > deno2)
               deno = deno1;
           else
               deno = deno2;
           end    
           m(i+3) =  m(i+2) - ((2*c)/deno);
           ae(i) = abs(((m(i+3) - m(i+2))/(m(i+3))));           
           if(abs(polyval(poly,m(i+3)))<=valf)
               flag = 1;
               break
           end
           if(ae(i)<= approxerror)
               flag = 2;
               break
           end    
       end
       fprintf('The value of the root is: %s\n',m(i+3));
       if flag == 2
          disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
       elseif flag==1
            disp('Stopping Criteria is convergence criteria for the function value');
       else
          disp('Stopping criteria is number of iterations');
       end
       y = @(x) poly(4)+poly(3)*x+poly(2)*x.^2+poly(1)*x.^3;
       fplot(y);
       grid on;
       plot(1:i,ae,'-o');
       grid on;
%       fplot(y);
%       figure;
%       plot(1:i,ae,'-o');       
    elseif(method ==2)
           rapproxerror = input('Enter stopping criteria\nRelative approximate error in r: ');
           sapproxerror = input('Enter stopping criteria\nRelative approximate error in s: ');
           itno = input('Maximum Iteration number: ');
           for k= 1:order/2
               r = [];
               s = [];
               r(1) = input('\nEnter r: ');
               s(1) = input('\nEnter s: ');           
               for j = 1:itno
                   q(1) = 1;
                   q(2) = -r(j);
                   q(3) = -s(j);
                   d(1) = poly(1);
                   d(2) = poly(2)+r(j)*d(1);
                   for i = 3:order+1
                       d(i) = poly(i) + r(j)*d(i-1) + s(j)*d(i-2);
                   end
                   e(1) = 0;
                   e(2) = d(1);
                   e(3) = d(2)+r(j)*e(2);
                   for i = 4:order+1
                       e(i) = d(i-1) + r(j)*e(i-1) + s(j)*e(i-2);
                   end
                   syms dr ds
                   eqn1 = e(order)*ds + e(order+1)*dr == -d(order+1);
                   eqn2 = e(order-1)*ds + e(order)*dr == -d(order);
                   sol = solve([eqn1, eqn2], [ds, dr]);
                   r(j+1) = r(j) + double(sol.dr);
                   s(j+1) = s(j) + double(sol.ds);
                   aer(j) = abs((r(j+1)-r(j))/r(j+1));
                   aes(j) = abs((s(j+1)-s(j))/s(j+1));
                   if(aer(j)<=rapproxerror && aes(j)<=sapproxerror)
                       flag = 2;
                       break;
                   end 
               end
               root1(1) = 0.5*(r(j)+sqrt(r(j)^2+(4*s(j))));
               root2(1) = 0.5*(r(j)-sqrt(r(j)^2+(4*s(j))));
               fprintf('The value of the root is: %f%+fi\n',real(root1(1)), imag(root1(1)));
               fprintf('The value of the root is: %f%+fi\n',real(root2(1)), imag(root2(1)));
               if flag == 2
                  disp('Stopping Criteria is convergence criterion for relative approximate errors in successive iterations');
               elseif flag==1
                    disp('Stopping Criteria is convergence criteria for the function value');
               else
                  disp('Stopping criteria is number of iterations'); 
               end
               y = @(x) poly(4)+poly(3)*x+poly(2)*x.^2+poly(1)*x.^3;
               plot(1:j,aer,'-o')
               grid on;
               figure;
               plot(1:j,aes,'-o')
               grid on;
               figure;
           end
           figure;
           fplot(y);
   end 
end