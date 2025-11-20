function [n, Fn, Ft] = getFunc( funcname )
switch funcname

    case 'DF1'
    n = 10;
    Fn = 2;
   
    H = @(t)0.75*sin(pi*t/2)+1.25;
    G = @(t)abs(sin(pi*t/2));
    f = @(x)x;
    g = @(x, t)1+sum((x-G(t)).^2);
    h = @(f,g,t)1-(f/g)^H(t);
    Fhandles = {@(x,t)f(x(1)), ...
                @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
    case 'DF2'
    n = 10;
    Fn = 2;
   
    G = @(t)abs(sin(pi*t/2));
    r= @(t)1+floor((n-1)*G(t));
    f = @(x)x;
    g = @(x, t)1+sum((x-G(t)).^2);
    Fhandles = {@(x,t)f(x(r(t))), ...
                @(x,t)(g(x(1:r(t)-1),t)+g(x(r(t)+1:n),t))*(1-sqrt(f(x(r(t)))/(g(x(1:r(t)-1),t)+g(x(r(t)+1:n),t))))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);   
    
    case 'DF3'
    n = 10;
    Fn = 2;
   
    G = @(t)sin(pi*t/2);
    H = @(t)G(t)+1.5;
    f = @(x)x;
    g = @(x, t)1+sum((x-G(t)-x(1)^H(t)).^2);
    Fhandles = {@(x,t)f(x(1)), ...
                @(x,t)g(x(2:end),t)*(1-(x(1)/g(x(2:n),t)).^H(t))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

    case 'DF4'
    n = 10;
    Fn = 2;
    
    a = @(t)sin(pi*t/2.0);
    b = @(t)1 + abs(cos(pi*t/2.0));
    H = @(t)a(t)+1.5; 
    g = @(x,t)1;
    for i=2:n
        g=@(x,t)g(x,t)+(x(i)-(a(t)*x(1)^2/i))^2;
    end  
    Fhandles = {@(x,t)g(x,t)*abs(x(1)-a(t))^H(t) , ...
                @(x,t)g(x,t)*abs(x(1)-a(t)-b(t))^H(t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
    
    case 'DF5'
    n = 10;
    Fn = 2;
    
    G = @(t)sin(pi*t/2);
    w =@(t)floor(10*G(t));
    g = @(x, t)1+sum((x-G(t)).^2);
    Fhandles = {@(x,t)g(x(2:n),t)*(x(1)+0.02*sin(w(t)*pi*x(1))) , ...
                @(x,t)g(x(2:n),t)*(1-x(1)+0.02*sin(w(t)*pi*x(1)))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
   
    
    case 'DF6'
    n = 10;
    Fn = 2;
    
    G = @(t)sin(pi*t/2);
    a = @(t)0.2+2.8*abs(G(t));
    y =@(x,t)x-G(t);
    g =@(x,t)1+ sum(abs(G(t))*y(x,t).^2-10*cos(2*pi*y(x,t))+10);
    Fhandles = {@(x,t)g(x(2:n),t)*(x(1)+0.1*sin(3*pi*x(1))).^a(t) , ...
                @(x,t)g(x(2:n),t)*(1-x(1)+0.1*sin(3*pi*x(1))).^a(t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
    
    case 'DF7'
    n = 10;
    Fn = 2;
    
    a = @(t)5*cos(pi*t/2);
    g =@(x,t)1+sum((x-1/(1+exp(a(t)*(x(1)-2.5)))).^2);
    Fhandles = {@(x,t)g(x(2:n),t)*(1+t)/x(1) , ...
                @(x,t)g(x(2:n),t)*x(1)/(1+t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);  
    
    case 'DF8'
    n = 10;
    Fn = 2;
    
    G = @(t)sin(pi*t/2);
    a = @(t)2.25+2*cos(2*pi*t);
    b = @(t)100*G(t).^2;
    g =@(x,t)1+ sum((x-G(t)*sin(power(4*pi*x(1),b(t)))/(1+abs(G(t)))).^2);
    Fhandles = {@(x,t)g(x(2:n),t)*(x(1)+0.1*sin(3*pi*x(1))) , ...
                @(x,t)g(x(2:n),t)*(1-x(1)+0.1*sin(3*pi*x(1))).^a(t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);   
    
   case 'DF9'
    n = 10;
    Fn = 2;
    
    N = @(t)1+floor(10*abs(sin(0.5*pi*t)));
    g =@(x,x1,t)1+ sum(x-cos(4*t+x(1)+x1)).^2;
    Fhandles = {@(x,t)g(x(2:n),x(1:n-1),t)*(x(1)+max(0,(0.1+0.5/N(t))*sin(2*N(t)*pi*x(1)))) , ...
                @(x,t)g(x(2:n),x(1:n-1),t)*(1-x(1)+max(0,(0.1+0.5/N(t))*sin(2*N(t)*pi*x(1)))) };
    Ft = @(x,t)cellfun(@(c)c(x,t),Fhandles);  
    
    case 'DF10'
    n = 10;
    Fn = 3; 
    
    H = @(t)2*cos(pi*t/2)+2.25;
    G = @(t)sin(pi*t/2);
    g = @(x,t)1+sum((x-sin(2*pi*(x(1)+x(2)))/(1+abs(G(t)))).^2);
    Fhandles = {@(x,t)g(x(3:n),t)*power(sin(0.5*pi*x(1)),H(t)) , ...
                @(x,t)g(x(3:n),t)*power(sin(0.5*pi*x(2)),H(t))*power(cos(0.5*pi*x(1)),H(t)), ...
                @(x,t)g(x(3:n),t)*power(cos(0.5*pi*x(2)),H(t))*power(cos(0.5*pi*x(1)),H(t))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);  
    
    
    case 'DF11'
    n = 10;
    Fn = 3; 
    
    G = @(t)abs(sin(pi*t/2));
    y= @(x,t)pi*G(t)/6.0+(pi/2-pi*G(t)/3.0)*x;
    g = @(x,t)1+G(t)+sum((x-0.5*G(t)*x(1)).^2);
    Fhandles = {@(x,t)g(x(3:n),t)*sin(y(x(1),t)), ...
                @(x,t)g(x(3:n),t)*sin(y(x(2),t))*cos(y(x(1),t)), ...
                @(x,t)g(x(3:n),t)*cos(y(x(2),t))*cos(y(x(1),t))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);  
    
    
   case 'DF12'
    n = 10;
    Fn = 3; 
    
    k=@(t)10*sin(pi*t);
    g=@(x,t)1+sum((x-sin(k(t)*x(1))).^2)+abs(sin(floor(k(t)*(2*x(1)-1))*pi/2)*sin(floor(k(t)*(2*x(2)-1))*pi/2));
    Fhandles = {@(x,t)g(x(3:n),t)*cos(0.5*pi*x(2))*cos(0.5*pi*x(1)), ...
                @(x,t)g(x(3:n),t)*sin(0.5*pi*x(2))*cos(0.5*pi*x(1)), ...
                @(x,t)g(x(3:n),t)*sin(0.5*pi*x(2))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);  
    
    case 'DF13'
    n = 10;
    Fn = 3; 
    
    G=@(t)sin(0.5*pi*t);
    p=@(t)floor(6*G(t));
    g=@(x,t)1+sum((x-G(t)).^2);
    Fhandles = {@(x,t)g(x(3:n),t)*cos(0.5*pi*x(1))^2, ...
                @(x,t)g(x(3:n),t)*cos(0.5*pi*x(2))^2, ...
                @(x,t)g(x(3:n),t)*sin(0.5*pi*x(1))^2+sin(0.5*pi*x(1))*cos(p(t)*pi*x(1))^2+sin(0.5*pi*x(2))^2+sin(0.5*pi*x(2))*cos(p(t)*pi*x(2))^2};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);    
    
    
    case 'DF14'
    n = 10;
    Fn = 3; 
    
    G=@(t)sin(0.5*pi*t);
    g=@(x,t)1+sum((x-G(t)).^2);
    y=@(x,t)0.5+G(t)*(x(1)-0.5);
    Fhandles = {@(x,t)g(x(3:n),t)*(1-y(x,t)+0.05*sin(6*pi*y(x,t))), ...
                @(x,t)g(x(3:n),t)*(1-x(2)+0.05*sin(6*pi*x(2)))*(y(x,t)+0.05*sin(6*pi*y(x,t))), ...
                @(x,t)g(x(3:n),t)*(x(2)+0.05*sin(6*pi*x(2)))*(y(x,t)+0.05*sin(6*pi*y(x,t)))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
    
    case 'FDA1'
    n = 10;
    Fn = 2;
   
    G = @(t)sin(pi*t/2);
    f = @(x)x;
    g = @(x, t)1+sum((x-G(t)).^2);
    Fhandles = {@(x,t)f(x(1)), ...
                @(x,t)1-sqrt(f(x(1)/g(x(2:end),t)))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
%     case 'FDA2'
%     n = 10;
%     Fn = 2;
%    
%     G = @(t)sin(pi*t/2);
%     H = @(t)0.7*sin(0.5*pi*t)+0.75;
%     f = @(x)x;
%     g = @(x, t)1+sum(x.^2);
%     h = @(x, t)sum((x-H(t)).^2);
%     Fhandles = {@(x,t)f(x(1)), ...
%                 @(x,t)g(x(2:n/2),t)*(1-f(x(1)/g(x(2:n/2),t))^(1/(H(t)+h(x(n/2+1:n)))))};
%     Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

    
    
    case 'FDA4'
    %% FDA4
    n = 10;
    Fn = 3;
    M = 3;

    G = @(t)abs(sin(pi*t/2));
    g = @(x,t)sum((x-G(t)).^2);
    Fhandles = {@(x,t)(1+g(x(M:n),t))*prod(cos(x(1:M-1)*pi/2)), ...
        @(x,t)(1+g(x(M:n),t))*prod(cos(x(1:M-2)*pi/2))*sin(x(M-1)*pi/2), ...
        @(x,t)(1+g(x(M:n),t))*sin(x(1)*pi/2)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    
case 'FDA5'
    %% FDA5
    n = 12;
    Fn = 3;
    M = 3;

    G = @(t)abs(sin(pi*t/2));
    g = @(x,t)G(t)+sum((x-G(t)).^2);
    F = @(t)1+100*sin(pi*t/2)^4;
    y = @(x,t)x.^F(t);
    Fhandles = {@(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-1), t)*pi/2)), ...
        @(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-2), t)*pi/2))*sin(y(x(M-1), t)*pi/2), ...
        @(x,t)(1+g(x(M:n),t))*sin(y(x(1),t)*pi/2)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'FDA5_iso'
    %% FDA5_iso
    n = 12;
    Fn = 3;
    M = 3;

    G = @(t)abs(sin(pi*t/2));
    B = .001;
    C = .05;
    g = @(x,t)G(t)+sum((min(0, floor(x-B)).*(G(t)*(B-x)/B-min(0,floor(C-x))*(1-G(t)).*(x-C)/(1-C))).^2);
    K = @(t)1+100*sin(pi*t/2)^4;
    y = @(x,t)x.^K(t);
    Fhandles = {@(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-1), t)*pi/2)), ...
        @(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-2), t)*pi/2))*sin(y(x(M-1), t)*pi/2), ...
        @(x,t)(1+g(x(M:n),t))*sin(y(x(1),t)*pi/2)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'FDA5_dec'
    %% FDA5_dec
    n = 12;
    Fn = 3;
    M = 3;

    G = @(t)abs(sin(pi*t/2));
    B = .001;
    C = .05;
    g = @(x,t)G(t)+sum(((floor(x-G(t)+B)*(1-C+(G(t)-B)/B)/(G(t)-B)+1/B+floor(G(t)+B-x)*(1-C+(1-G(t)-B)/B)/(1-G(t)-B)).*(abs(x-G(t))-B)+1).^2);
    K = @(t)1+100*sin(pi*t/2)^4;
    y = @(x,t)x.^K(t);
    Fhandles = {@(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-1), t)*pi/2)), ...
        @(x,t)(1+g(x(M:n),t))*prod(cos(y(x(1:M-2), t)*pi/2))*sin(y(x(M-1), t)*pi/2), ...
        @(x,t)(1+g(x(M:n),t))*sin(y(x(1),t)*pi/2)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'DIMP2'
    %% DIMP2   <--- [-2,2]
    n = 10;
    Fn = 2;


    G = @(i,t)sin(pi*t/2+2*pi*(i/(n+1))).^2;
    f = @(x)x;
    g = @(x, t)1+2*(n-1)+sum((4*x-2-G(2:n,t)).^2-2*cos(3*pi*(4*x-2-G(2:n,t))));
    h = @(f,g)1-sqrt(f/g);
    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t))};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'dMOP2'
    %% dMOP2
    n = 10;
    Fn = 2;


    H = @(t)0.75*sin(pi*t/2)+1.25;
    G = @(t)sin(pi*t/2);
    f = @(x)x;
    g = @(x, t)1+9*sum((x-G(t)).^2);
    h = @(f,g,t)1-(f/g)^H(t);
    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'dMOP2_iso'
    %% dMOP2_iso
    n = 10;
    Fn = 2;

    
    H = @(t)0.75*sin(pi*t/2)+1.25;
    B = .001;
    C = .05;
    G = @(t)sin(pi*t/2);
    f = @(x)x;
    g = @(x,t)1+9*sum((min(0, floor(x-B)).*(G(t)*(B-x)/B-min(0,floor(C-x))*(1-G(t)).*(x-C)/(1-C))).^2);
    h = @(f,g,t)1-(f/g)^H(t);
    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'dMOP2_dec'
    %% dMOP2_dec
    n = 10;
    Fn = 2;

    
    H = @(t)0.75*sin(pi*t/2)+1.25;
    B = .001;
    C = .05;
    G = @(t)sin(pi*t/2);
    f = @(x)x;
    g = @(x,t)1+9*sum((((floor(x-G(t)+B)*(1-C+(G(t)-B)/B)/(G(t)-B)+1/B+floor(G(t)+B-x)*(1-C+(1-G(t)-B)/B)/(1-G(t)-B)).*(abs(x-G(t))-B)+1)-G(t)).^2);
    h = @(f,g,t)1-(f/g)^H(t);
    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'dMOP3'
    %% dMOP3
    n = 10;
    Fn = 2;

    
    H = @(t)0.75*sin(pi*t/2)+1.25;
    G = @(t)sin(pi*t/2);
    f = @(x)x;
    g = @(x, t)1+9*sum((x-G(t)).^2);
    h = @(f,g,t)1-sqrt(f/g);
    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'HE2'
    %% HE2
    n = 30;
    Fn = 2;
    
    

    H = @(t)0.75*sin(pi*t/2)+1.25;
    f = @(x)x;
    g = @(x, t)1+9*sum(x)/(n-1);
    h = @(f,g,t)1-sqrt(f/g)^H(t)-(f/g)^H(t)*sin(10*pi*f);

    Fhandles = {@(x,t)f(x(1)), ...
        @(x,t)g(x(2:end),t)*h(f(x(1)), g(x(2:end), t), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
    F = @(x)Ft(x,t);

case 'HE7'
    %% HE7 <---- [-1,1]
    n = 10;
    Fn = 2;


    H = @(t)0.75*sin(pi*t/2)+1.25;
    f = @(x)x(1)+2/numel(3:2:n)*sum((2*x(3:2:n)-1-(0.3*x(1)^2*cos(24*pi*x(1)+4*(3:2:n)*pi/n)+0.6*x(1)).*cos(6*pi*x(1)+(3:2:n)*pi/n)).^2);
    g = @(x)2-sqrt(x(1))+2/numel(2:2:n)*sum((2*x(2:2:n)-1-(0.3*x(1)^2*cos(24*pi*x(1)+4*(2:2:n)*pi/n)+0.6*x(1)).*sin(6*pi*x(1)+(2:2:n)*pi/n)).^2);
    h = @(f,g,t)1-(f/g)^H(t);

    Fhandles = {@(x,t)f(x), ...
        @(x,t)g(x)*h(f(x), g(x), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);

case 'HE9'
    %% HE9
    n = 10;
    Fn = 2;

    H = @(t)0.75*sin(pi*t/2)+1.25;
    f = @(x)x(1)+2/numel(3:2:n)*sum((x(3:2:n)-sin(6*pi*x(1)+(3:2:n)*pi/n)).^2);
    g = @(x)2-x(1)^2+2/numel(2:2:n)*sum((x(2:2:n)-sin(6*pi*x(1)+(2:2:n)*pi/n)).^2);
    h = @(f,g,t)1-(f/g)^H(t);

    Fhandles = {@(x,t)f(x), ...
        @(x,t)g(x)*h(f(x), g(x), t)};
    Ft = @(x,t)cellfun(@(c)c(x,t), Fhandles);
otherwise
    disp([funcname ' Not Found!'])
end