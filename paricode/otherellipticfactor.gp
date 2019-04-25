\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ compare mahler measure of nonreciprocal families and L-values of elliptic factors
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

read("mahler.gp");
read("mahlergeneral.gp");
read("jGenus1.gp");
read("comparej.gp");


Alist=[1,x^2+x+1,\
x^2+x+1,\
(x^2+x+1)^2,\
1,1,(x^2+x+1),(x^2+x+1),\
(x^2+x+1)^2,(x^2+x+1)^2,\
(x^2+x+1)^2,(x^2+x+1)^2\
];
Blist=[x^3+k*x^2+k*x+1,x^3+k*x^2+k*x+1,\
k*x^2+k*x,\
k*x^3+k*x+(2*k+2)*x^2,\
x^4+1+k*x^3+k*x+(2*k-4)*x^2,x^4+1+k*x^3+k*x+2*k*x^2,\
x^4+1+k*x^3+k*x+(2*k-4)*x^2,x^4+1+k*x^3+k*x+2*k*x^2,\
x^4+1+k*x^3+k*x+(2*k-4)*x^2,x^4+1+k*x^3+k*x+2*k*x^2,\
2*x^4+2+k*x^3+k*x+(2*k-6)*x^2,2*x^4+2+k*x^3+k*x+(2*k-2)*x^2\
];
Clist=[x^3,x*(x^2+x+1),\
x*(x^2+x+1),\
(x^2+x+1)^2,\
x^4,x^4,(x^2+x+1)*x^2,(x^2+x+1)*x^2,\
(x^2+x+1)^2,(x^2+x+1)^2,\
(x^2+x+1)^2,(x^2+x+1)^2\
];

ellfactorlist1=[y^2 + (2*x^2 + (1 - k)*x + 1)*y + x*(x-1)*(x+1)^2,\
y^3 + (-x + 3)*y^2 - (x^2 + (-k + 1)*x - 3)*y + (x+1)*(x^2+x+1),\
x^4 - x^3 + (k*y - 1)*x + (-y^2 + 1),\
0,\
y^2 + (x^2 + (k - 2)*x + 1)*y + x*(x - 1)^2,\
(x + 1)*y^2 + (x^2 + (2-k)*x + 1)*y + (x^2 + x),\
y^2 + (k - 2)*x*y - (x-1)^2*(x^2+x+1),\
(x^2 + x + 1)*y^2 + (x^2 + (-k + 2)*x + 1)*y + (x^2 + x + 1),\
0,\
0,\
0,\
0];
ellfactorlist2=[y^2 + ((k - 3)*x + 1)*y + x^3,\
y^2 + (x^2 + (-k + 3)*x - 1)*y + x^2*(x - 1),\
x^2 + (-y^2 + k*y + 2)*x + 1,\
x^4 + 4*x^3 + 6*x^2 + (-k*y + 4)*x + (y^2 + 2*y + 1),\
y^2 + (4 - k)*x*y + x*(x^2 - 2*x + 1),\
y^2 + (x^2 + (k - 4)*x + 1)*y + x^2,\
y^2 + ((k - 4)*x + 2)*y - (x-1)*(x+1)^2,\
(x + 1)*y^2 + (x^2 + (-k + 4)*x + 1)*y + (x^2+x),\
(y^2 + y + 1)*x^2 + (y^2 + (k - 4)*y + 1)*x + (y^2 + y + 1),\
x^4 + 4*x^3 + (y + 6)*x^2 + ((-k + 4)*y + 4)*x + (y^2 + 2*y + 1),\
y^2 + (k - 8)*x*y + (x^2+1)^2,\
(x+1)^2*y^2 + (2*x^2 + (-k + 8)*x + 2)*y + (x+1)^2];

{
calcjInvariantOther() =
my(vec=vector(#Alist));
for(num=1,#Alist,
my(A=Alist[num]);
my(C=Clist[num]);
my(B=Blist[num]);
my(cl=calcEllipticFactor(B^2-4*A*C));
my(E = ellinit(cl));
vec[num]=[E.j,A*y^2+B*y+C];
);
return(vec);
}

\\ write results to the tex file
filepath="compareother.tex";


{
\\ find elliptic curve defined by tempered polynomial
\\ which is isomorphic to the elliptic factors in the Jacobian of the genus 2 curve
compare() =
my(jlist1=calcjInvariantOther());
my(jlist2=getj(1));
for(i=1,#jlist1,
  for(j=1,#jlist2,
    my(ji=subst(jlist1[i][1],k,k));
    my(jj=subst(jlist2[j][1],k,k));
    degni=poldegree(numerator(ji));
    degnj=poldegree(numerator(jj));
    degdi=poldegree(denominator(ji));
    degdj=poldegree(denominator(jj));
    if(degni==degnj && degdi==degdj,
      compfunc=compareRationalFunction(ji,jj);
      if(compfunc!=0,print(i);print(j);print(compfunc);
        \\k2=a*k1+b
        my(a=-polcoeff(compfunc,1,k1)/polcoeff(compfunc,1,k2));
        my(b=-subst(subst(compfunc,k1,0),k2,0)/polcoeff(compfunc,1,k2));
        write(filepath,"==================================");
        write(filepath,i);
        write(filepath,j);
        my(bstr="");
        if(b<0,bstr=b);
        if(b>0,bstr=Str("+",b));
        my(kStr);
        if(a==-1,kStr=Str("-k",bstr),kStr=Str("k",bstr));
        \\write(filepath,kStr);
        write(filepath,polylist[j]+(a*k+b)*x*y);
        \\comparemahlerother(i,jlist2[j][2],a,b);
      );
    );
  );
);
}


{
\\ compare mahler measure between P(x-1,y) and the elliptic factors
comparemahlerother(i,P2,a,b) =
for(kk=-40,40,
  my(kkk=kk/2);
  my(P1=Alist[i]*y^2+Blist[i]*y+Clist[i]);
  my(delta=subst(Blist[i],k,kkk)^2-4*Alist[i]*Clist[i]);
  cl=calcEllipticFactor2(delta);
  a2=eval(cl[2]);
  a4=eval(cl[4]);
  a6=eval(cl[5]);
  disc = eval(18*a2*a4*a6 - 4*a2^3*a6 + a2^2*a4^2 - 4*a4^3 - 27*a6^2);
  P=subst(subst(P1,k,kkk),x,x-1);
  my(d=gcd(polcoeff(P,2,y),polcoeff(P,1,y)));
  if(d!=1,P=polcoeff(P,2,y)/d*y^2+polcoeff(P,1,y)/d*y+polcoeff(P,0,y)/d);
  print(P);
  my(m1=calcmahlergeneral());
  print(m1);
  P=subst(P2,k,a*kkk+b);
  print(P);
  \\ reciprocal polynomial
  if(checkreciprocal(P)!=0,
    my(d=gcd(polcoeff(P,2,y),polcoeff(P,1,y)));
    if(d!=1,P=polcoeff(P,2,y)/d*y^2+polcoeff(P,1,y)/d*y+polcoeff(P,0,y)/d);
  );
  my(m2=calcmahlergeneral());
  print(m2);
  write(filepath,Str("k = ",kkk));
  if(m2!=0,write(filepath,m1/m2));
  if(disc!=0 && (1.0*kkk)%1==0,
    my(E = ellinit(eval(cl)));
    my(s=lfun(E, 0, 1)/m1);
    write(filepath,Str("s = ",s));
  );
);
}


{
calcnumofroots(i)=
for(kk=-30,30,
  my(poly=subst(Blist[i],k,kk)^2-4*Alist[i]*Clist[i]);
  if(poly%(x+1)^2==0,poly=poly/(x+1)^2);
  if(poly!=0,
    my(rouc=calcrouc(poly));
    print(#rouc);
  );
);
}

{
\\ roots on unit circle
calcrouc(poly) =
my(rod=polroots(poly));
my(rouc=listcreate(length(rod)));
my(j=1);
\\ find all the roots on the unit circle exp(2*Pi*I*x)
for(i=1,length(rod),
  if(abs(abs(rod[i])-1)<10^(-20),
  listput(rouc,imag(log(rod[i]))/2/Pi,j);
  j=j+1)
);
listsort(rouc);
return(rouc);
}


{
numofroots()=
for(i=1,#Alist,
  delta=(Blist[i]^2-4*Alist[i]*Clist[i]);
  if(delta%(x+1)^2==0,delta=delta/(x+1)^2);
  print("======================================");
  for(kk=-40,40,
    if(simplify(subst(delta,k,kk))==0,return);
    my(rlist=polroots(subst(delta,k,kk)));
    my(num=0);
    for(i=1,#rlist,
      if(abs(rlist[i]+1)<1,num=num+1);
    );
    if(num%2==0, print(kk));
  );
);
}


{
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ find linear relation between the Mahler measures of P(x-1,y)
\\ and the polynomials in ellfactorlist1 and ellfactorlist2
\\ the results:
\\ 1  m1=m3 k<=0  2*m1=m2+m3 k>=7
\\ 2  m1=m3 k<=2  2*m1=m2+m3 k>=8
\\ 5  m1=m3 k<=-8 2*m1=m2+m3 k>=9
\\ 6  m1=m3 k<=-1 2*m1=m2+m3 k>=17
\\ 7  m1=m3 k<=-7 2*m1=m2+m3 k>=8
\\ 8  m1=m3 k<=1  2*m1=m2+m3 k>=16
\\ 11 m1=m3 k<=-4
\\ 12 m1=m3 k<=4
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
comparethreemahler()=
for(i=1,#Alist,
  write(filepath,Str("============",i,"============"));
  my(P1=Alist[i]*y^2+Blist[i]*y+Clist[i]);
  my(P2=ellfactorlist1[i]);
  my(P3=ellfactorlist2[i]);
  if(i==11||i==12,P2=P1);
  for(j=-50,50,
    jj=j;
    P=subst(subst(P1,k,jj),x,x-1);
    P=simplifyP(P);
    my(m1=calcmahlergeneral());
    my(m2=1);
    if(P2!=0,
      P=subst(P2,k,jj);
      \\ reciprocal polynomial
      if(checkreciprocal(P)!=0,
        P=simplifyP(P);
      );
      m2=calcmahlergeneral();
    );

    P=subst(P3,k,jj);
    \\ reciprocal polynomial
    if(checkreciprocal(P)!=0,
      P=simplifyP(P);
    );
    my(m3=calcmahlergeneral());
    print(m3/m1);
    dep=lindep([m1,m2,m3],20);
    if(dep[1]<5 && dep[2]<5 && dep[3]<5,
      print(j);
      print(dep);
      write(filepath,jj);
      write(filepath,dep);
    );
  );
);
}

{
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ find linear relation between the Mahler measures of P(x-1,y) and L-values of elliptic factors
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
comparemahlerlvalue()=
write(filepath,"\\documentclass{amsart}");
write(filepath,"\\usepackage{longtable}");
write(filepath,"\\usepackage{amssymb}");
write(filepath,"\\begin{document}");
for(i=1,#Alist,
  write(filepath,"\\begin{longtable}{llllll}");
  write(filepath,"\\hline");
  write(filepath,"$k$ & $c_1$ & $c_2$ & $c_3$ & $N_{E}$ & $N_{F}$\\\\");
  write(filepath,"\\hline");
  my(P1=Alist[i]*y^2+Blist[i]*y+Clist[i]);
  for(j=-50,50,
    jj=j;
    P=subst(subst(P1,k,jj),x,x-1);
    P=simplifyP(P);
    my(m1=calcmahlergeneral());

    delta=Blist[i]^2-4*Alist[i]*Clist[i];
    cl1=subst(calcEllipticFactor(delta),k,jj);
    a2=eval(cl1[2]);
    a4=eval(cl1[4]);
    a6=eval(cl1[5]);
    disc1 = eval(18*a2*a4*a6 - 4*a2^3*a6 + a2^2*a4^2 - 4*a4^3 - 27*a6^2);

    cl2=subst(calcEllipticFactor2(delta),k,jj);
    a2=eval(cl2[2]);
    a4=eval(cl2[4]);
    a6=eval(cl2[5]);
    disc2 = eval(18*a2*a4*a6 - 4*a2^3*a6 + a2^2*a4^2 - 4*a4^3 - 27*a6^2);
    my(l1=1);
    my(l2=1);
    my(N1=0);
    my(N2=0);
    if(disc1 !=0,
      my(E1 = ellinit(eval(cl1)));
      l1=lfun(E1, 0, 1);
      my(egr1=ellglobalred(E1));
      N1=egr1[1];
    );
    if(disc2 != 0,
      my(E2 = ellinit(eval(cl2)));
      l2=lfun(E2, 0, 1);
      my(egr2=ellglobalred(E2));
      N2=egr2[1];
    );
    dep=lindep([m1,l1,l2],20);
    if(dep[1]!=0 && dep[3]!=0 && denominator(dep[1]/dep[3])<20,
      write(filepath, Str(jj, " & ", dep[1], " & ", dep[2], " & ", dep[3], " & ", N1, " & ", N2,  "\\\\"));
    );
  );
  write(filepath,"\\hline");
  write(filepath,"\\end{longtable}");
);
write(filepath,"\\end{document}");
}


{
simplifyP(P)=
  my(d=gcd(polcoeff(P,2,y),polcoeff(P,1,y)));
  if(d!=1,P=polcoeff(P,2,y)/d*y^2+polcoeff(P,1,y)/d*y+polcoeff(P,0,y)/d);
  return(P);
}
