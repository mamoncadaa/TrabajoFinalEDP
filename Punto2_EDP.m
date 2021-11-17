%Limites
inf= 0;
sup = 1;
%Epsilon
e = 0.11;
errori=ones(6,1);
hs=ones(6,1);
alpha=ones(6,1);
for j=1:6
    hC = @(j) 2.^(-j);
    h = hC(j);
    X = (inf:h:sup);
    n = size(X,2)-1;
    Su = @(x) 1 + x + (exp(x/e)-1)./(exp(1/e)-1);    
    upa =  1 - (h/(2*e));
    ua = -2;
    ujs =  1 + (h/(2*e));
    fj = - h^2/e;
    C = fj;
    C(1) = fj -  (1 + (h/(2*e)));
    C(end) =  fj -  3*(1 - (h/(2*e)));
    Uaprx = ones(n+1,1);
    Uaprx(1) = 1;
    Uaprx(end) = 3;
    Uaprx([2:end-1]) = tridiag(ua,upa,ujs,C);
    Uexact = Su(X)';
    Error = Uaprx-Uexact;
    
    Tabla = [X' Uaprx Uexact Error];
    Nvariable = {'x', 'U aprox', 'U exacta', 'Error absoluto'};
    Tabla = array2table(Tabla, 'VariableNames',Nvariable);
    
    
    cont=1;
    integrala=0;
    for i=inf:h:sup-h
       interp=@(x)abs((((Uaprx(cont)-Uaprx(cont+1))/(i-(i+h)))*(x-i)+Uaprx(cont))-(1 + x + (exp(x/e)-1)./(exp(1/e)-1))).^2;
       integrala=integral(interp,i,i+h)+integrala  ;
       cont=cont+1;
    end
    errori(j)=integrala^(1/2);
    hs(j)=h;

    tablalatex.data = [X' Uaprx Uexact Error];
    tablalatex.tableColLabels = Nvariable
    latex = latexTable(tablalatex);
   
    plot(X', Uaprx)
    hold on
    plot(X', Uexact)
    
    legend('h=1/2', ...
    ' h=1/4', ...
    ' h=1/8', ...
    ' h=1/16',...
    ' h=1/32',...
    ' h=1/64',...
    ' Exacta')
    title("epsilon =0.0001" )
end
alpha(1)=nan
for i=2:6;
    alpha(i)=log(errori(i)/errori(i-1))/log(hs(i)/hs(i-1));
end

Tabla2 = [hs errori alpha];
Nvariable2 = {'h', 'e_h', 'alpha_h'};
Tabla2 = array2table(Tabla2, 'VariableNames',Nvariable2);
disp(Tabla2)
tablalatex.data = [hs errori alpha];
tablalatex.tableColLabels = Nvariable2
latex = latexTable(tablalatex);   
%loglog(hs,errori)         
%xlabel('h')
%ylabel('E_h')     
    


    
    
   
